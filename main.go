package main

import (
	"bufio"
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

// To do
func alignClusters(refFilePath string, clusterBamPath string) {
	f, _ := os.Open(clusterBamPath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	// read ref file

	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}
		// extract corresponding ref segment
		fmt.Printf("%s\n", rec.Name)
		// do stuff with the read & ref

	}
}

var segmentSize int

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	svTag = sam.NewTag("SV")

	var signalingReads map[string][]IndexPair

	segmentSize = findAverageSegmentSize(*bamFile, 1000, 1000000)
	fmt.Printf("Segment size = %d\n", segmentSize)
	svStore, ciStore := readVcfFiltered(*vcfFile, "del")
	fmt.Printf("\nDone reading vcf %d %d\n", len(svStore.svMap), len(ciStore.ciList))
	/*
		signalingReads = matchAlignmentsToSVs(*bamFile, ciStore)
		WriteMapToFile(signalingReads, "signalingReads40.txt")
	*/
	signalingReads = restoreSignalingReads("signalingReads40.txt")
	clusters := constructClusterFile(*bamFile, signalingReads, svStore, ciStore, "unmapped")
	for key, value := range clusters {
		fmt.Printf("I found %d reads for interval %d\n", len(value), key)
	}

	// To do
	alignClusters(*refFile, "unmapped_cluster.bam")

}
