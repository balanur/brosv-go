package main

import (
	"fmt"
	"io"
	"log"
	"os"
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

	l, r = getRefParts(currentCI, currentSV.Type)
	refL := ref.getChr(currentSV.Chromosome).content[l.head : l.tail+1]
	refR := ref.getChr(currentSV.Chromosome).content[r.head : r.tail+1]
	read := string(rec.Seq.Expand())

	result := align(len(read), refL, refR, read, currentSV.Type)

	if (result.identityL >= 0.95 && result.identityR >= 0.95) || (result.identityL == -1 && result.identityR >= 0.95) || (result.identityL >= 0.95 && result.identityR == -1) {
		/*
			cigarL, _ := computeCIGAR(result.aL, result.bL, result.cL)
			cigarR, rlen := computeCIGAR(result.aR, result.bR, result.cR)

			var cigar sam.Cigar
			svlen := (r.head + result.rbp) - (l.head + result.lbp) - 1
			if currentSV.Type == "INV" {
				svlen -= len(read) - rlen
			}
			if svlen > 0 {
				svreg := sam.NewCigarOp(sam.CigarDeletion, svlen)
				cigar = append(cigarL, svreg)
				cigar = append(cigar, cigarR...)
			} else {
				cigar = append(cigarL, cigarR...)
			}
		*/
		rec.Pos = result.pos + l.head
		rec.Ref.SetName(currentSV.Chromosome)
		//rec.Cigar = cigar
		newAux, _ := sam.NewAux(lbpTag, result.lbp+l.head)
		rec.AuxFields = append(rec.AuxFields, newAux)
		newAux, _ = sam.NewAux(rbpTag, result.rbp+r.head)
		rec.AuxFields = append(rec.AuxFields, newAux)
		return rec, true
	}
	return nil, false
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
	return true
}

func getRefParts(currentCI Interval, svtype string) (Interval, Interval) {
	var l, r Interval
	if currentCI.side == leftCI {
		l = currentCI
		r = ciStore.ciList[rightCIs[currentCI.svId]]
	} else {
		l = ciStore.ciList[leftCIs[currentCI.svId]]
		r = currentCI
	}
	if svtype == "INV" {
		if currentCI.side == leftCI {
			r.head -= segmentSize
		} else {
			l.tail -= segmentSize
		}
	}
	return l, r
}
