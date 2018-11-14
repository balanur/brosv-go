package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

func alignClusters(refFilePath string, clusterBamPath string, outfilePath string, svStore SVStore, ciStore CIStore) {
	f, _ := os.Open(clusterBamPath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	g1, _ := os.Create(outfilePath + ".bam")
	defer g1.Close()

	alWriter, _ := bam.NewWriter(g1, bamReader.Header(), 0)
	defer alWriter.Close()

	ref := readReference(refFilePath)
	var alignments []*sam.Record

	var wg sync.WaitGroup
	wg.Add(*threads)
	channels := make([]chan *sam.Record, *threads)
	var resultLock sync.Mutex
	for threadIndex := 0; threadIndex < *threads; threadIndex++ {
		channels[threadIndex] = make(chan *sam.Record, 2000)
		go func(tIndex int) {
			defer wg.Done()
			for rec := range channels[tIndex] {
				alignedRec, flag := alignSingleRead(svStore, ciStore, ref, rec)
				if flag {
					resultLock.Lock()
					alignments = append(alignments, alignedRec)
					resultLock.Unlock()
				}
				fmt.Printf("Aligned %d reads\r", len(alignments))
			}
		}(threadIndex)
	}

	readIndex := 0
	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}

		channels[readIndex%(*threads)] <- rec
		readIndex++
		if readIndex%100 == 0 {
			fmt.Printf("Distributed %d, at %s %d\t\t\r", readIndex, rec.Ref.Name(), rec.Pos)
		}
	}

	for i := 0; i < *threads; i++ {
		close(channels[i])
	}

	fmt.Println("\nWaiting on threads")
	wg.Wait()

	// write alignments to bam
	fmt.Printf("Writing alignment results to file\n")
	for i := 0; i < len(alignments); i++ {
		alWriter.Write(alignments[i])
	}
}

func alignSingleRead(svStore SVStore, ciStore CIStore, ref Genome, rec *sam.Record) (*sam.Record, bool) {
	var l, r Interval
	ciIndex := auxValue(rec.AuxFields.Get(svTag))
	currentCI := ciStore.ciList[ciIndex]
	currentSV := svStore.get(currentCI.svId)

	if currentCI.side == leftCI {
		l = currentCI
		r = ciStore.ciList[rightCIs[currentCI.svId]]
	} else {
		l = ciStore.ciList[leftCIs[currentCI.svId]]
		r = currentCI
	}
	refL := ref.getChr(currentSV.Chromosome).content[l.head : l.tail+1]
	refR := ref.getChr(currentSV.Chromosome).content[r.head : r.tail+1]
	read := string(rec.Seq.Expand())

	var result AlignmentResult
	result = align(len(read), refL, refR, read)
	if rec.Flags&sam.Reverse != 0 {
		rec.Flags = sam.Paired + sam.ProperPair + sam.Reverse
	} else {
		rec.Flags = sam.Paired + sam.ProperPair + sam.MateReverse
	}

	if (result.identityL >= 0.95 && result.identityR >= 0.95) || (result.identityL == -1 && result.identityR >= 0.95) || (result.identityL >= 0.95 && result.identityR == -1) {
		cigarL, _ := computeCIGAR(result.aL, result.bL, result.cL)
		cigarR, _ := computeCIGAR(result.aR, result.bR, result.cR)

		var cigar sam.Cigar
		dellen := (r.head + result.rbp) - (l.head + result.lbp) - 1
		if dellen > 0 {
			delreg := sam.NewCigarOp(sam.CigarDeletion, dellen)
			cigar = append(cigarL, delreg)
			cigar = append(cigar, cigarR...)
		} else {
			cigar = append(cigarL, cigarR...)
		}
		rec.Pos = result.pos + l.head
		rec.Ref.SetName(currentSV.Chromosome)
		rec.Cigar = cigar
		newAux, _ := sam.NewAux(lbpTag, result.lbp+l.head)
		rec.AuxFields = append(rec.AuxFields, newAux)
		newAux, _ = sam.NewAux(rbpTag, result.rbp+r.head)
		rec.AuxFields = append(rec.AuxFields, newAux)
		return rec, true
	}
	return nil, false
}

func extractBreakpointResults(bamFilePath string, outfile string, ciStore CIStore, svStore SVStore) {
	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %v:\n", err)
	}

	g, _ := os.Create(outfile)
	defer g.Close()

	writer := bufio.NewWriter(g)

	supportedSVs := make(map[string]int)
	lbp := make(map[string][]int)
	rbp := make(map[string][]int)
	count := 0
	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}

		ciIndex := auxValue(rec.AuxFields.Get(svTag))
		currentCI := ciStore.ciList[ciIndex]

		if _, exist := supportedSVs[currentCI.svId]; exist {
			supportedSVs[currentCI.svId]++
			lbp[currentCI.svId] = append(lbp[currentCI.svId], auxValue(rec.AuxFields.Get(lbpTag)))
			rbp[currentCI.svId] = append(rbp[currentCI.svId], auxValue(rec.AuxFields.Get(rbpTag)))
		} else {
			supportedSVs[currentCI.svId] = 1
			lbp[currentCI.svId] = append(lbp[currentCI.svId], auxValue(rec.AuxFields.Get(lbpTag)))
			rbp[currentCI.svId] = append(rbp[currentCI.svId], auxValue(rec.AuxFields.Get(rbpTag)))
			count++
		}
	}

	writer.WriteString("Number of supported SVs: " + strconv.Itoa(count) + "\n")
	fmt.Printf("Number of supported SVs: %d\n", count)
	for k, v := range supportedSVs {
		writer.WriteString("SVID: " + k + " Chromosome: " + svStore.svMap[k].Chromosome + " numofReads: " + strconv.Itoa(v) + "\n")
		left := ciStore.ciList[leftCIs[k]]
		right := ciStore.ciList[rightCIs[k]]
		writer.WriteString("Left: " + strconv.Itoa(left.head) + ", " + strconv.Itoa(left.tail))
		writer.WriteString(" Right: " + strconv.Itoa(right.head) + ", " + strconv.Itoa(right.tail) + "\n\n")
		for i := range lbp[k] {
			writer.WriteString(strconv.Itoa(lbp[k][i]) + "\t" + strconv.Itoa(rbp[k][i]) + "\n")
		}
		writer.WriteString("\n\n")
	}
	writer.Flush()
}

func isValidSplit(cigar []sam.CigarOp) bool {
	max := 0
	index := 0
	// find main del
	for i := 0; i < len(cigar); i++ {
		if cigar[i].Type() == sam.CigarDeletion && max < cigar[i].Len() {
			max = cigar[i].Len()
			index = i
		}
	}
	mlen := 0
	// check left
	for i := 0; i < index; i++ {
		if cigar[i].Type() == sam.CigarMatch {
			mlen += cigar[i].Len()
		}
	}
	if mlen < 5 {
		return false
	}
	mlen2 := 0
	// check right
	for i := index + 1; i < len(cigar); i++ {
		if cigar[i].Type() == sam.CigarMatch {
			mlen2 += cigar[i].Len()
		}
	}
	if mlen2 < 5 {
		return false
	}
	if mlen+mlen < 50 {
		return false
	}
	return true
}

func restoreRefs(filePath string) map[string]*sam.Reference {
	seenRefs := make(map[string]*sam.Reference)
	f, _ := os.Open(filePath)
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		id, _ := strconv.Atoi(words[1])
		len, _ := strconv.Atoi(words[2])
		ref, _ := sam.NewReference(words[0], "", "", len, nil, nil)
		ref.id = id
		seenRefs[words[0]] = ref
	}

	return seenRefs
}

func alignContigs(refFilePath string, clusterBamPath string, contigfile string, outfilePath string, ciStore CIStore, svStore SVStore, seenRefs map[string]*sam.Reference) {
	f, _ := os.Open(contigfile)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	f2, _ := os.Open(clusterBamPath)
	defer f.Close()
	bamReader, _ := bam.NewReader(f2, *threads)
	defer bamReader.Close()

	g, _ := os.Create(outfilePath + ".bam")
	//g, _ := os.Create("temp/test.bam")
	defer g.Close()
	alWriter, _ := bam.NewWriter(g, bamReader.Header(), 0)
	defer alWriter.Close()

	ref := readReference(refFilePath)
	var refL, refR string
	var currentCI, l, r Interval
	var currentSV SV
	var ciIndex int
	count := 0
	allcount := 0

	for scanner.Scan() {
		line := scanner.Text()
		// processing new CI
		if line[0] == 'c' {
			words := strings.Fields(scanner.Text())
			ciIndex, _ = strconv.Atoi(words[1])
			fmt.Printf("Processing ci %d reads mapped so far %d / %d\n", ciIndex, count, allcount)
			currentCI = ciStore.ciList[ciIndex]
			//currentSV = svStore.get(currentCI.svId)
			if currentCI.side == leftCI {
				l = currentCI
				r = ciStore.ciList[rightCIs[currentCI.svId]]
			} else {
				l = ciStore.ciList[leftCIs[currentCI.svId]]
				r = currentCI
			}
			refL = ref.getChr(currentSV.Chromosome).content[l.head : l.tail+1]
			refR = ref.getChr(currentSV.Chromosome).content[r.head : r.tail+1]
		} else if line[0] == '>' {
			continue
		} else {
			allcount++
			contig := scanner.Text()
			var result AlignmentResult
			var flags sam.Flags
			// Reverse compelement the left side
			if currentCI.side == leftCI {
				contig = Complement(Reverse(contig))
				result = align(len(contig), refL, refR, contig)
				flags = sam.Reverse
			} else {
				result = align(len(contig), refL, refR, contig)
				flags = sam.MateReverse
			}

			if (result.identityL >= 0.95 && result.identityR >= 0.95) || (result.identityL == -1 && result.identityR >= 0.95) || (result.identityL >= 0.95 && result.identityR == -1) {
				cigarL, _ := computeCIGAR(result.aL, result.bL, result.cL)
				cigarR, _ := computeCIGAR(result.aR, result.bR, result.cR)

				var cigar sam.Cigar
				dellen := (r.head + result.rbp) - (l.head + result.lbp) - 1
				if dellen > 0 {
					delreg := sam.NewCigarOp(sam.CigarDeletion, dellen)
					cigar = append(cigarL, delreg)
					cigar = append(cigar, cigarR...)
				} else {
					cigar = append(cigarL, cigarR...)
				}
				chr := ref.getChr(currentSV.Chromosome)
				ref, _ := sam.NewReference(chr.title, "", "", len(chr.content), nil, nil)
				ref.SetID(seenRefs[ref.Name()].ID())

				qual := make([]byte, len(contig))
				for i := range qual {
					qual[i] = 34
				}
				rec, err := sam.NewRecord(strconv.Itoa(count), ref, nil, result.pos+l.head, -1, len(contig), 0, cigar, []byte(contig), qual, []sam.Aux{})
				if err != nil {
					log.Fatalf("problem %v\n", err)
				}
				rec.Flags = flags
				// lbp & rbp tags
				newAux, _ := sam.NewAux(lbpTag, result.lbp+l.head)
				rec.AuxFields = append(rec.AuxFields, newAux)
				newAux2, _ := sam.NewAux(rbpTag, result.rbp+r.head)
				rec.AuxFields = append(rec.AuxFields, newAux2)
				// sv tag
				newAux3, _ := sam.NewAux(svTag, ciIndex)
				rec.AuxFields = append(rec.AuxFields, newAux3)
				err = alWriter.Write(rec)
				count++
			}
		}
	}
	fmt.Printf("All contigs vs. mapped ones: %d / %d\n", allcount, count)
}
