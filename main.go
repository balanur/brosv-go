package main

import (
	"bufio"
	"bytes"
	"flag"
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
	"github.com/brentp/vcfgo"
)

var (
	mode = flag.Int("mode", 0, "Running mode.\n"+
		"0: Generate Signaling Reads\n"+
		"1: Generate Custers\n")
	vcfFile = flag.String("vcf", "", "vcf input file")
	bamFile = flag.String("bam", "", "bam input file")
	refFile = flag.String("ref", "", "reference file")
	sr      = flag.String("sr", "", "Signaling Reads file")
	threads = flag.Int("threads", 0, "number of threads to use (0 = auto)")
	help    = flag.Bool("help", false, "display help")
)

var svTag sam.Tag
var segmentSize int
var leftCIs, rightCIs map[string]int

func readVcf(fileName string) (SVStore, CIStore) {
	return readVcfFiltered(fileName, "")
}

func readVcfFiltered(fileName string, filter string) (SVStore, CIStore) {
	svStore := NewSVStore()
	ciStore := NewCIStore()

	f, _ := os.Open(fileName)
	rdr, err := vcfgo.NewReader(f, false)

	if err != nil {
		panic(err)
	}

	for {
		variant := rdr.Read()
		if variant == nil {
			break
		}
		if !strings.Contains(variant.Id(), filter) {
			continue
		}
		var tempSV SV
		tempSV.Start = int(variant.Pos)
		endPosition, _ := variant.Info().Get("END")
		tempSV.End = endPosition.(int)
		tempSV.Chromosome = variant.Chromosome
		svType, _ := variant.Info().Get("SVTYPE")
		tempSV.Type = svType.(string)
		tempSV.id = variant.Id()

		svStore.add(tempSV)

		var leftInterval Interval
		if ciPos, err := variant.Info().Get("CIPOS"); err == nil {
			leftInterval.head = tempSV.Start + ciPos.([]int)[0]
			leftInterval.tail = tempSV.Start + ciPos.([]int)[1]
			leftInterval.svId = tempSV.id
			leftInterval.side = leftCI
		}

		var rightInterval Interval
		if ciEnd, err := variant.Info().Get("CIEND"); err == nil {
			rightInterval.head = tempSV.End + ciEnd.([]int)[0]
			rightInterval.tail = tempSV.End + ciEnd.([]int)[1]
			rightInterval.svId = tempSV.id
			rightInterval.side = rightCI
		}

		if leftInterval.head != leftInterval.tail {
			leftInterval.head -= segmentSize
			ciStore.add(svStore, leftInterval)
		}

		if rightInterval.head != rightInterval.tail {
			rightInterval.tail += segmentSize
			ciStore.add(svStore, rightInterval)
		}
	}
	fmt.Printf("SV is %d / %d", len(ciStore.ciList), len(svStore.svMap)*2)
	return svStore, ciStore
}

func isSignaling(record *sam.Record) bool {
	flags := record.Flags
	pos := record.Pos
	matePos := record.MatePos

	if flags&sam.Paired == 0 {
		return false
	}

	if flags&sam.Unmapped != 0 {
		return false
	}

	// Mate is unmapped.
	if flags&sam.MateUnmapped != 0 {
		return true
	}

	// Mate is in another chromosome
	if record.Ref.Name() != record.MateRef.Name() {
		return true
	}

	// Same direction with mate
	if flags&sam.Reverse != 0 && flags&sam.MateReverse != 0 { // --
		return true
	}
	if flags&sam.Reverse == 0 && flags&sam.MateReverse == 0 { // ++
		return true
	}

	// Read placed before its mate -+
	if flags&sam.Reverse != 0 && flags&sam.MateReverse == 0 && pos <= matePos {
		return true
	}

	// Read placed after its mate -+
	if flags&sam.Reverse == 0 && flags&sam.MateReverse != 0 && pos > matePos {
		return true
	}

	return false
}

func matchAlignmentsToSVs(bamFilePath string, ciStore CIStore) map[string][]IndexPair {
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
				if isSignaling(rec) {
					intersectingIntervals := ciStore.findIntersectingIntervals(rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())
					for _, intervalIndex := range intersectingIntervals {
						pair := IndexPair{mappingOri: getMappingOri(rec), pairNumber: getPairNumber(rec), ciIndex: intervalIndex}
						//fmt.Printf("Signaling Read for ci: %d, %d %d rec.name: %s\n", intervalIndex, getMappingOri(rec), getPairNumber(rec), rec.Name)
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

func findAverageSegmentSize(bamFilePath string, thr int, N int) int {
	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	i := 0
	sum := 0
	for {
		rec, err := bamReader.Read()
		if err == io.EOF || i >= N {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}

		if rec.Flags&sam.Paired != 0 && rec.Flags&sam.ProperPair != 0 && rec.TempLen > 0 && rec.TempLen < thr {
			sum += rec.TempLen
			i++
		}
	}

	return sum / i
}

func restoreSignalingReads(filePath string) map[string][]IndexPair {
	signalingReads := make(map[string][]IndexPair)

	f, _ := os.Open(filePath)
	defer f.Close()

	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		var indices []IndexPair
		for i := 1; i < len(words); i += 3 {
			indexCI, _ := strconv.Atoi(words[i])
			pairNumber, _ := strconv.Atoi(words[i+1])
			mappingOri, _ := strconv.Atoi(words[i+2])
			indices = append(indices, IndexPair{ciIndex: indexCI, pairNumber: pairNumber, mappingOri: mappingOri})
		}
		signalingReads[words[0]] = indices
	}
	return signalingReads
}

func constructClusterFile(bamFilePath string, signalingReads map[string][]IndexPair, svStore SVStore, ciStore CIStore, outputBamFilePath string) map[int][]*sam.Record {
	// First go over bamfile, take every pair, put them into a map where key is indexCI, value is array of reads.
	ciReadMap := make(map[int][]*sam.Record)

	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	g1, _ := os.Create(outputBamFilePath + "_cluster.bam")
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
				if indices[i].pairNumber != getPairNumber(rec) {
					interval := ciStore.get(indices[i].ciIndex)
					/*sv := svStore.get(interval.svId)
					fmt.Printf("I'm mate: %d/%d SV Id: %s Chr: %s\n", interval.head, interval.tail, sv.id, sv.Chromosome)
					fmt.Scanf("\n")*/
					if interval.side == leftCI && indices[i].mappingOri == 2 {
						continue
					}
					if interval.side == rightCI && indices[i].mappingOri == 1 {
						continue
					}

					newAux, _ := sam.NewAux(svTag, indices[i].ciIndex)
					rec.AuxFields = append(rec.AuxFields, newAux)
					unmappedWriter.Write(rec)
					ciReadMap[indices[i].ciIndex] = append(ciReadMap[indices[i].ciIndex], rec)
				}
			}
		}

		readIndex++
		if readIndex%1000000 == 0 {
			fmt.Printf("Reads at %s %d %d\r", rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())
		}
	}

	return ciReadMap
}

func alignClusters(refFilePath string, clusterBamPath string, outfilePath string, svStore SVStore, ciStore CIStore) {
	f, _ := os.Open(clusterBamPath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	f2, _ := os.Create(outfilePath)
	defer f2.Close()

	writer := bufio.NewWriter(f2)

	g1, _ := os.Create("alignment.bam")
	defer g1.Close()

	alWriter, _ := bam.NewWriter(g1, bamReader.Header(), 0)
	defer alWriter.Close()

	ref := readReference(refFilePath)
	var refL, refR string
	var currentCI, l, r Interval
	ciIndex := -1
	var side, before int
	count := 0
	allcount := 0

	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}
		temp := auxValue(rec.AuxFields.Get(svTag))
		if ciIndex != temp {
			ciIndex = temp
			fmt.Printf("Processing ci %d\n", ciIndex)
			currentCI = ciStore.ciList[ciIndex]
			currentSV := svStore.get(currentCI.svId)
			if currentCI.side == leftCI {
				side = 1
				l = currentCI
				r = ciStore.ciList[rightCIs[currentCI.svId]]
			} else {
				side = 2
				l = ciStore.ciList[leftCIs[currentCI.svId]]
				r = currentCI
			}
			refL = ref.getChr(currentSV.Chromosome).content[l.head : l.tail+1]
			refR = ref.getChr(currentSV.Chromosome).content[r.head : r.tail+1]
		}
		read := string(rec.Seq.Expand())
		if rec.Flags&sam.Reverse == 0 && currentCI.side == leftCI {
			read = Complement(Reverse(read))
		}
		if rec.Flags&sam.Reverse != 0 && currentCI.side == rightCI {
			read = Complement(Reverse(read))
		}
		//var lbp, rbp, acc int
		//var al string
		var result AlignmentResult
		initAligner()
		if currentCI.side == leftCI {
			//lbp, rbp, acc, al = alignExtend(refL, refR, read)
			result = align(refL, refR, read)
			rec.Flags = sam.Paired + sam.ProperPair + sam.Reverse
		} else {
			//rbp, lbp, acc, al = alignExtend(Reverse(refR), Reverse(refL), Reverse(read))
			//rbp = len(refR) - rbp - 1
			//lbp = len(refL) - lbp - 1
			result = align(refL, refR, read)
			rec.Flags = sam.Paired + sam.ProperPair + sam.MateReverse
		}

		if result.identityL+result.identityR >= 1.9 {
			if before != ciIndex {
				writer.WriteString(strconv.Itoa(ciIndex) + " " + strconv.Itoa(side) + ":\n")
				before = ciIndex
			}

			writer.WriteString("lbp: " + strconv.Itoa(result.lbp+l.head) + " rbp: " + strconv.Itoa(result.rbp+r.head) + "\n")
			writer.WriteString(result.aL + "\n" + result.cL + "\n" + result.bL + "\nidentityL: " + strconv.FormatFloat(result.identityL, 'f', 3, 64) + "\n\n")
			writer.WriteString(result.aR + "\n" + result.cR + "\n" + result.bR + "\nidentityR: " + strconv.FormatFloat(result.identityR, 'f', 3, 64) + "\n\n")
			cigarL, len1 := computeCIGAR(result.aL, result.bL, result.cL)
			cigarR, len2 := computeCIGAR(result.aR, result.bR, result.cR)

			var cigar sam.Cigar
			dellen := (r.head + result.rbp) - (l.head + result.lbp) - 1
			if dellen > 0 {
				delreg := sam.NewCigarOp(sam.CigarDeletion, dellen)
				cigar = append(cigarL, delreg)
				cigar = append(cigar, cigarR...)
			} else {
				cigar = append(cigarL, cigarR...)
			}
			rec.Pos = result.lbp + l.head
			rec.Cigar = cigar
			alWriter.Write(rec)
			if len1+len2 != len(read) {
				fmt.Printf("%s\n", read)
				fmt.Printf("%s\n", cigar.String())
				fmt.Printf("%s\n%s\n%s\n\n", result.aL, result.cL, result.bL)
				fmt.Printf("%s\n%s\n%s\n\n", result.aR, result.cR, result.bR)
				fmt.Printf("%d %d %d\n", len1+len2, len(cigar.String()), len(read))
				break
			}
			count++
		}
		allcount++
	}
	writer.Flush()
	fmt.Printf("%d out of %d reads are mapped\n", count, allcount)
}

func alignExtend(refL string, refR string, read string) (int, int, int, string) {
	var lbp, rbp int
	var a, b, c, al bytes.Buffer
	i := 0
	match := 0
	mismatch := 0

	//left- 10 exact match
	kmer := read[0:10]
	loc := strings.Index(refL, kmer)
	//extend left
	if loc != -1 {
		a.WriteString(read[0:10])
		b.WriteString("||||||||||")
		c.WriteString(refL[loc : loc+10])
		i = i + 10     // read iterator
		lbp = loc + 10 // ref iterator
		match = 10
		for mismatch <= 1 && i < len(read) && lbp < len(refL) {
			if read[i] != refL[lbp] {
				mismatch++
				b.WriteString(" ")
			} else {
				b.WriteString("|")
				match++
			}
			a.WriteString(read[i : i+1])
			c.WriteString(refL[lbp : lbp+1])
			i++
			lbp++
		}
	} else { // no match at left
		al.WriteString(a.String())
		al.WriteString("\n")
		al.WriteString(b.String())
		al.WriteString("\n")
		al.WriteString(c.String())
		al.WriteString("\n")
		return -1, -1, 0, al.String()
	}
	a.WriteString("xxxxxxxxxx")
	c.WriteString("xxxxxxxxxx")
	b.WriteString("          ")
	//right- 10 exact match
	if i+10 > len(read) {
		kmer = read[i:len(read)]
	} else {
		kmer = read[i : i+10]
	}
	loc = strings.Index(refR, kmer)
	//extend right
	if loc != -1 {
		a.WriteString(kmer)
		for j := 0; j < len(kmer); j++ {
			b.WriteString("|")
		}
		c.WriteString(refR[loc : loc+len(kmer)])
		i = i + len(kmer)
		rbp = loc + len(kmer)
		mismatch = 0
		match += len(kmer)
		for mismatch <= 1 && i < len(read) && rbp < len(refR) {
			if read[i] != refR[rbp] {
				mismatch++
				b.WriteString(" ")
			} else {
				b.WriteString("|")
				match++
			}
			a.WriteString(read[i : i+1])
			c.WriteString(refR[rbp : rbp+1])
			i++
			rbp++
		}
	} else {
		rbp = -1
	}

	acc := 100 * match / len(read)
	al.WriteString(a.String())
	al.WriteString("\n")
	al.WriteString(b.String())
	al.WriteString("\n")
	al.WriteString(c.String())
	al.WriteString("\n")
	return lbp, rbp, acc, al.String()
}
func linkLeftRightCIs(svStore SVStore, ciStore CIStore) {
	leftCIs = make(map[string]int)
	rightCIs = make(map[string]int)

	for i, interval := range ciStore.ciList {
		if interval.side == leftCI {
			leftCIs[interval.svId] = i
		} else {
			rightCIs[interval.svId] = i
		}
	}
}

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	svTag = sam.NewTag("SV")

	//var signalingReads map[string][]IndexPair

	segmentSize = findAverageSegmentSize(*bamFile, 1000, 1000000)
	fmt.Printf("Segment size = %d\n", segmentSize)
	svStore, ciStore := readVcfFiltered(*vcfFile, "del")
	linkLeftRightCIs(svStore, ciStore)

	fmt.Printf("\nDone reading vcf %d %d\n", len(svStore.svMap), len(ciStore.ciList))
	/*
		signalingReads = matchAlignmentsToSVs(*bamFile, ciStore)
		WriteMapToFile(signalingReads, "signalingReads40.txt")
		signalingReads = restoreSignalingReads("signalingReads40.txt")
		clusters := constructClusterFile(*bamFile, signalingReads, svStore, ciStore, "unmapped")
		for key, value := range clusters {
			fmt.Printf("I found %d reads for interval %d\n", len(value), key)
		}
	*/
	alignClusters(*refFile, "unmapped_cluster_sorted.bam", "trash.txt", svStore, ciStore)

}
