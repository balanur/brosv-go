package main

import (
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

func isSignaling(record *sam.Record, svType SVType) bool {
	flags := record.Flags
	pos := record.Pos
	matePos := record.MatePos
	cigar := record.Cigar.String()

	if flags&sam.Paired == 0 {
		return false
	}

	if flags&sam.Unmapped != 0 {
		return false
	}

	// Mate is unmapped
	if flags&sam.MateUnmapped != 0 {
		return false
	}

	// Mate is in another chromosome
	if record.Ref.Name() != record.MateRef.Name() {
		return false
	}

	// Read placed before its mate -+ //tandup
	if flags&sam.Reverse != 0 && flags&sam.MateReverse == 0 && pos <= matePos {
		return false
	}

	// Read placed after its mate -+ //tandup
	if flags&sam.Reverse == 0 && flags&sam.MateReverse != 0 && pos > matePos {
		return false
	}

	// Clipped Alignments Pattern
	r, _ := regexp.Compile("[1-9][0-9]S")

	if svType == inv {

		// Same direction with mate
		if flags&sam.Reverse != 0 && flags&sam.MateReverse != 0 { // --
			return true
		}
		if flags&sam.Reverse == 0 && flags&sam.MateReverse == 0 { // ++
			return true
		}

		// just split in inv region
		if r.MatchString(cigar) {
			return true
		}
	}

	if svType == del {
		// Insert size (pairs mapped too closer or farther than expected)
		max := float64(segmentSize + 3*variance)
		if math.Abs(float64(pos-matePos)) > max {
			if r.MatchString(cigar) {
				return true
			}
		}
		// Just clipped (not mapping too farther than expected)
		if r.MatchString(cigar) {
			return true
		}

		if strings.Index(cigar, "H") != -1 {
			return true
		}
	}
	return false
}

func extractSignalingInCI(bamFilePath string, outputBamFilePath string, ciStore CIStore, svType SVType) {

	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	g, _ := os.Create(outputBamFilePath)
	defer g.Close()

	bamWriter, _ := bam.NewWriter(g, bamReader.Header(), 0)
	defer bamWriter.Close()

	var wg sync.WaitGroup
	wg.Add(*threads)
	channels := make([]chan *sam.Record, *threads)
	var writeLock sync.Mutex

	for threadIndex := 0; threadIndex < *threads; threadIndex++ {
		channels[threadIndex] = make(chan *sam.Record, 2000)
		go func(tIndex int) {
			defer wg.Done()

			for rec := range channels[tIndex] {

				if isSignaling(rec, svType) {

					intersectingIntervals := ciStore.findIntersectingIntervals(rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())

					for _, intervalIndex := range intersectingIntervals {
						newAux, _ := sam.NewAux(svTag, intervalIndex)
						rec.AuxFields = append(rec.AuxFields, newAux)
						writeLock.Lock()
						bamWriter.Write(rec)
						writeLock.Unlock()
						rec.AuxFields = rec.AuxFields[:len(rec.AuxFields)-1]
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

}

func setBreakpointTags(bamFilePath string, outputBamFilePath string, ciStore CIStore) {
	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	g1, _ := os.Create(outputBamFilePath + "cluster_withbp.bam")
	defer g1.Close()

	bamWriter, _ := bam.NewWriter(g1, bamReader.Header(), 0)
	defer bamWriter.Close()

	readIndex := 0
	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}

		ciIndex := auxValue(rec.AuxFields.Get(svTag))
		currentCI := ciStore.ciList[ciIndex]
		cigar := rec.Cigar.String()

		m := strings.Index(cigar, "M")
		s := strings.Index(cigar, "S")
		h := strings.Index(cigar, "H")
		var val int

		if m < s || m < h {
			m, _ = strconv.Atoi(cigar[0:m])
			val = m + rec.Pos
		} else {
			if s != -1 {
				m, _ = strconv.Atoi(cigar[s+1 : m])
			} else if h != -1 {
				m, _ = strconv.Atoi(cigar[h+1 : m])
			}
			val = rec.Pos
		}

		// eliminate insignificant splits
		if m >= 10 {
			if currentCI.side == leftCI {
				newAux, _ := sam.NewAux(lbpTag, val)
				rec.AuxFields = append(rec.AuxFields, newAux)
			} else {
				newAux, _ := sam.NewAux(rbpTag, val)
				rec.AuxFields = append(rec.AuxFields, newAux)
			}

			bamWriter.Write(rec)
			rec.AuxFields = rec.AuxFields[:len(rec.AuxFields)-1]
		}
		readIndex++
		if readIndex%1000000 == 0 {
			fmt.Printf("Reads at %s %d %d\r", rec.Ref.Name(), rec.Pos, rec.Pos+rec.Len())
		}
	}
}
