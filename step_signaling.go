package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

func isSignaling(record *sam.Record) (SVType, int) {
	flags := record.Flags
	pos := record.Pos
	matePos := record.MatePos
	cigar := record.Cigar.String()

	if flags&sam.Paired == 0 {
		return none, 0
	}

	if flags&sam.Unmapped != 0 {
		return none, 0
	}

	// Mate is unmapped
	if flags&sam.MateUnmapped != 0 {
		return none, 0
	}

	// Mate is in another chromosome
	if record.Ref.Name() != record.MateRef.Name() {
		//return true, false - not for dels
		return none, 0
	}

	// Same direction with mate
	if flags&sam.Reverse != 0 && flags&sam.MateReverse != 0 { // --
		return inv, 0
	}
	if flags&sam.Reverse == 0 && flags&sam.MateReverse == 0 { // ++
		return inv, 0
	}

	// Read placed before its mate -+
	if flags&sam.Reverse != 0 && flags&sam.MateReverse == 0 && pos <= matePos {
		return tandup, 0
	}

	// Read placed after its mate -+
	if flags&sam.Reverse == 0 && flags&sam.MateReverse != 0 && pos > matePos {
		return tandup, 0
	}

	// Discard hard clipped alignments
	if strings.Index(cigar, "H") != -1 {
		return none, 0
	}

	// Clipped Alignments
	r, _ := regexp.Compile("[1-9][0-9]S")

	// Insert size (pairs mapped too closer or farther than expected)
	max := float64(segmentSize + 3*variance)
	if math.Abs(float64(pos-matePos)) > max {
		if r.MatchString(cigar) {
			return del, 2
		} else {
			return del, 1
		}
	}

	// Just clipped (not mapping too farther than expected)
	if r.MatchString(cigar) {
		return del, 2
	}
	return none, 0
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

				svType, whichCi := isSignaling(rec)

				if svType == del {
					intersectingIntervals := ciStore.findIntersectingIntervals(rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())

					// write for the interval itself
					for _, intervalIndex := range intersectingIntervals {
						pair := IndexPair{mappingOri: getMappingOri(rec), pairNumber: getPairNumber(rec), ciIndex: intervalIndex}
						mapLock.Lock()
						signalingReads[rec.Name] = append(signalingReads[rec.Name], pair)
						mapLock.Unlock()
					}
					// write for other CI of SV as well
					if whichCi == 2 {
						for _, intervalIndex := range intersectingIntervals {
							currentCI := ciStore.ciList[intervalIndex]
							var otherCI int
							if currentCI.side == leftCI {
								otherCI = rightCIs[currentCI.svId]
							} else {
								otherCI = leftCIs[currentCI.svId]
							}

							pair := IndexPair{mappingOri: getMappingOri(rec), pairNumber: getPairNumber(rec), ciIndex: otherCI}
							mapLock.Lock()
							signalingReads[rec.Name] = append(signalingReads[rec.Name], pair)
							mapLock.Unlock()
						}
					}
				} else if svType == inv { // to be continued...
					intersectingIntervals := ciStore.findIntersectingIntervals(rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())
					for _, intervalIndex := range intersectingIntervals {
						pair := IndexPair{mappingOri: getMappingOri(rec), pairNumber: getPairNumber(rec), ciIndex: intervalIndex}
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
			writer.WriteString(" " + strconv.Itoa(v[i].ciIndex) + " " + strconv.Itoa(v[i].pairNumber) + " " + strconv.Itoa(v[i].mappingOri))
		}
		writer.WriteString("\n")
	}

	writer.Flush()
}
