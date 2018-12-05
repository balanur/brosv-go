// THIS IS A BACKUP FILE IN CASE I NEED THESE OLD FUNCTIONS LATER

package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

// --------------------------- Organizer functions ---------------------------
func clusteringMode() {
	fmt.Printf("Running in mode 1 - Clustering\n")
	fmt.Printf("Clustering for dels < 10,000\n")
	signalingReads := restoreSignalingReads(*workdir + *sr)
	constructClusterFile(*bamFile, signalingReads, svStore, ciStore, *workdir)
	setBreakpointTags(path.Join(*workdir, "cluster.bam"), *workdir, ciStore)
}

func assemblyMode() {
	fmt.Printf("Running in mode 2 - Assembling clusters\n")
	assembleReads(path.Join(*workdir, "cluster_sorted.bam"))
}

func alignmentMode() {
	fmt.Printf("Running in mode 3 - Aligning clusters\n")
	//seenRefs := restoreRefs(path.Join(*workdir, "seenRefs.txt"))
	//alignContigs(*refFile, path.Join(*workdir, "cluster.bam"), path.Join(*workdir, "simu.contigs"), path.Join(*workdir, "alignment_contig"), ciStore, svStore, seenRefs)
	alignClusters(*refFile, path.Join(*workdir, "cluster_sorted.bam"), path.Join(*workdir, "alignment40"), svStore, ciStore)
	extractBreakpointResults(path.Join(*workdir, "alignment40.bam"), path.Join(*workdir, "supportedSVs.txt"), ciStore, svStore)
}

// --------------------------- Signaling Read Extraction ---------------------------
func matchAlignmentsToSVs(bamFilePath string, ciStore CIStore, svType SVType) map[string][]IndexPair {
	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	signalingReads := make(map[string][]IndexPair)

	var wg sync.WaitGroup
	wg.Add(*threads)
	channels := make([]chan *sam.Record, *threads)
	var mapLock sync.Mutex
	for threadIndex := 0; threadIndex < *threads; threadIndex++ {
		channels[threadIndex] = make(chan *sam.Record, 2000)
		go func(tIndex int) {
			defer wg.Done()

			for rec := range channels[tIndex] {

				if isSignaling(rec, svType) {

					intersectingIntervals := ciStore.findIntersectingIntervals(rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())

					for _, intervalIndex := range intersectingIntervals {
						pair := IndexPair{mappingOri: getMappingOri(rec), pairNumber: getPairNumber(rec), ciIndex: intervalIndex, secondary: getSecondary(rec)}
						mapLock.Lock()
						signalingReads[rec.Name] = append(signalingReads[rec.Name], pair)
						mapLock.Unlock()
					}
				}
			}
		}(threadIndex)
	}

	fmt.Println("Distributing to threads")

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
		if readIndex%100000 == 0 {
			fmt.Printf("Distributed %d, at %s %d\t\t\r", readIndex, rec.Ref.Name(), rec.Pos)
		}
	}

	for i := 0; i < *threads; i++ {
		close(channels[i])
	}
	fmt.Println("Waiting on threads")
	wg.Wait()

	return signalingReads
}

func WriteMapToFile(signalingReads map[string][]IndexPair, filePath string) {
	f, _ := os.Create(filePath)
	defer f.Close()

	writer := bufio.NewWriter(f)

	for k, v := range signalingReads {
		writer.WriteString(k)
		// If value type is integer array
		for i := range v {
			writer.WriteString(" " + strconv.Itoa(v[i].ciIndex) + " " + strconv.Itoa(v[i].pairNumber) + " " + strconv.Itoa(v[i].mappingOri) + " " + strconv.Itoa(v[i].secondary))
		}
		writer.WriteString("\n")
	}

	writer.Flush()
}

// --------------------------- Clustering ---------------------------

func restoreSignalingReads(filePath string) map[string][]IndexPair {
	signalingReads := make(map[string][]IndexPair)

	f, _ := os.Open(filePath)
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		var indices []IndexPair
		for i := 1; i < len(words); i += 4 {
			indexCI, _ := strconv.Atoi(words[i])
			pairNumber, _ := strconv.Atoi(words[i+1])
			mappingOri, _ := strconv.Atoi(words[i+2])
			secondary, _ := strconv.Atoi(words[i+3])
			indices = append(indices, IndexPair{ciIndex: indexCI, pairNumber: pairNumber, mappingOri: mappingOri, secondary: secondary})
		}
		signalingReads[words[0]] = indices
	}
	return signalingReads
}

func constructClusterFile(bamFilePath string, signalingReads map[string][]IndexPair, svStore SVStore, ciStore CIStore, outputBamFilePath string) (map[int][]*sam.Record, map[string]*sam.Reference) {
	// First go over bamfile, take every pair, put them into a map where key is indexCI, value is array of reads.
	ciReadMap := make(map[int][]*sam.Record)
	seenReferences := make(map[string]*sam.Reference)

	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	g1, _ := os.Create(outputBamFilePath + "cluster.bam")
	defer g1.Close()

	unmappedWriter, _ := bam.NewWriter(g1, bamReader.Header(), 0)
	defer unmappedWriter.Close()

	readIndex := 0
	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}

		if indices, ok := signalingReads[rec.Name]; ok {
			for i := range indices {
				if indices[i].pairNumber == getPairNumber(rec) && indices[i].secondary == getSecondary(rec) {
					newAux, _ := sam.NewAux(svTag, indices[i].ciIndex)
					rec.AuxFields = append(rec.AuxFields, newAux)
					unmappedWriter.Write(rec)
					rec.AuxFields = rec.AuxFields[:len(rec.AuxFields)-1]

					ciReadMap[indices[i].ciIndex] = append(ciReadMap[indices[i].ciIndex], rec)
					seenReferences[rec.Ref.Name()] = rec.Ref
				}
			}
		}

		readIndex++
		if readIndex%1000000 == 0 {
			fmt.Printf("Reads at %s %d %d\r", rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())
		}
	}

	return ciReadMap, seenReferences
}

func WriteClusterToFile(clusters map[int][]*sam.Record, filePath string, ciStore CIStore, svStore SVStore) {
	f, _ := os.Create(filePath)
	defer f.Close()
	writer := bufio.NewWriter(f)
	writer.WriteString("Number of clusters(CIs) having reads:" + strconv.Itoa(len(clusters)) + "\n")
	for k, v := range clusters {
		ci := ciStore.ciList[k]
		chr := svStore.svMap[ci.svId].Chromosome
		writer.WriteString(chr + " " + strconv.Itoa(k) + " " + strconv.Itoa(len(v)) + "\n")
	}
	writer.Flush()
}

func WriteRefsToFile(seenRefs map[string]*sam.Reference, filePath string) {
	f, _ := os.Create(filePath)
	defer f.Close()
	writer := bufio.NewWriter(f)

	for k, v := range seenRefs {
		writer.WriteString(k + " " + strconv.Itoa(v.ID()) + " " + strconv.Itoa(v.Len()) + "\n")
	}
	writer.Flush()
}

func findSVnum(clusterfile string, ciStore CIStore) {
	f, _ := os.Open(clusterfile)
	defer f.Close()
	existSV := make(map[string]int)
	count := 0

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		id, _ := strconv.Atoi(words[1])
		currentCI := ciStore.ciList[id]
		if _, exist := existSV[currentCI.svId]; exist {
			existSV[currentCI.svId] = 2
		} else {
			count++
			existSV[currentCI.svId] = 1
		}
	}
	fmt.Printf("number of SVs having reads %d\n", len(existSV))
	fmt.Printf("number of SVs having reads %d\n", count)
}

// --------------------------- Assembly ---------------------------

func assembleReads(clusterBamPath string) {
	f, _ := os.Open(clusterBamPath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	g, _ := os.Create("tempfile.fastq")
	defer g.Close()
	writer := bufio.NewWriter(g)

	g2, _ := os.Create(*workdir + "simu.contigs")
	defer g2.Close()
	writer2 := bufio.NewWriter(g2)

	currentCI := -1
	for {
		rec, err := bamReader.Read()
		var ci int
		if err != io.EOF {
			ci = auxValue(rec.AuxFields.Get(svTag))
		}

		if currentCI == -1 {
			currentCI = ci
		}
		if currentCI != ci || err == io.EOF {

			writer.Flush()
			// do the assembly
			cmd := exec.Command("rm", "-rf", "velv")
			cmd.Run()
			cmd = exec.Command("velveth", "velv/", "60", "-fastq", "-short", "tempfile.fastq", "-strand_specific")
			cmd.Run()
			cmd = exec.Command("velvetg", "velv/", "-unused_reads", "yes", "-min_contig_lgth", "100", "-cov_cutoff", "1")
			cmd.Run()

			os.Remove("tempfile.fastq")
			g, _ = os.Create("tempfile.fastq")
			defer g.Close()
			writer = bufio.NewWriter(g)

			// write results
			contigs, _ := os.Open("velv/contigs.fa")
			defer contigs.Close()
			scanner := bufio.NewScanner(contigs)
			writer2.WriteString("ci " + strconv.Itoa(currentCI))
			for scanner.Scan() {
				line := scanner.Text()
				if line[0] == '>' {
					writer2.WriteString("\n")
					writer2.WriteString(line + "\n")
				} else {
					writer2.WriteString(line)
				}
			}
			writer2.WriteString("\n")

			if err == io.EOF {
				break
			}
			currentCI = ci
		}
		// add read sequence to fasta (read name, seq, ori, qual)
		writer.WriteString("@" + rec.Name + "\n")
		seq := string(rec.Seq.Expand())
		if rec.Flags&sam.Reverse != 0 {
			seq = Reverse(Complement(seq))
		}
		writer.WriteString(seq + "\n")
		writer.WriteString("+\n")
		qual := string(formatQual(rec.Qual))
		writer.WriteString(qual + "\n")
	}
	writer2.Flush()
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
		// Processing new CI
		if line[0] == 'c' {
			words := strings.Fields(scanner.Text())
			ciIndex, _ = strconv.Atoi(words[1])

			fmt.Printf("Processing ci %d reads mapped so far %d / %d\n", ciIndex, count, allcount)

			currentCI = ciStore.ciList[ciIndex]
			currentSV = svStore.get(currentCI.svId)

			l, r = getRefParts(currentCI, currentSV.Type)
			refL = ref.getChr(currentSV.Chromosome).content[l.head : l.tail+1]
			refR = ref.getChr(currentSV.Chromosome).content[r.head : r.tail+1]
		} else if line[0] == '>' {
			continue
		} else {
			allcount++
			contig := scanner.Text()
			var result AlignmentResult

			result = align(len(contig), refL, refR, contig, currentSV.Type)

			if (result.identityL >= 0.95 && result.identityR >= 0.95) || (result.identityL == -1 && result.identityR >= 0.95) || (result.identityL >= 0.95 && result.identityR == -1) {
				cigarL, _ := computeCIGAR(result.aL, result.bL, result.cL)
				cigarR, _ := computeCIGAR(result.aR, result.bR, result.cR)

				var cigar sam.Cigar
				var op sam.CigarOpType
				svlen := (r.head + result.rbp) - (l.head + result.lbp) - 1
				if svlen > 0 {
					if currentSV.Type == "del" {
						op = sam.CigarDeletion
					}
					svreg := sam.NewCigarOp(op, svlen)
					cigar = append(cigarL, svreg)
					cigar = append(cigar, cigarR...)
				} else {
					cigar = append(cigarL, cigarR...)
				}
				chr := ref.getChr(currentSV.Chromosome)
				ref, _ := sam.NewReference(chr.title, "", "", len(chr.content), nil, nil)
				//ref.SetID(seenRefs[ref.Name()].ID())

				qual := make([]byte, len(contig))
				for i := range qual {
					qual[i] = 34
				}
				rec, _ := sam.NewRecord(strconv.Itoa(count), ref, nil, result.pos+l.head, -1, len(contig), 0, cigar, []byte(contig), qual, []sam.Aux{})
				// lbp & rbp tags
				newAux, _ := sam.NewAux(lbpTag, result.lbp+l.head)
				rec.AuxFields = append(rec.AuxFields, newAux)
				newAux2, _ := sam.NewAux(rbpTag, result.rbp+r.head)
				rec.AuxFields = append(rec.AuxFields, newAux2)
				// sv tag
				newAux3, _ := sam.NewAux(svTag, ciIndex)
				rec.AuxFields = append(rec.AuxFields, newAux3)
				alWriter.Write(rec)
				count++
			}
		}
	}
	fmt.Printf("All contigs vs. mapped ones: %d / %d\n", allcount, count)
}

func restoreRefs(filePath string) map[string]*sam.Reference {
	seenRefs := make(map[string]*sam.Reference)
	f, _ := os.Open(filePath)
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		//id, _ := strconv.Atoi(words[1])
		len, _ := strconv.Atoi(words[2])
		ref, _ := sam.NewReference(words[0], "", "", len, nil, nil)
		//ref.SetID(id)
		seenRefs[words[0]] = ref
	}

	return seenRefs
}

// --------------------------- Alignment ---------------------------

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
