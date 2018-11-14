package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

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

		// Skip hard clipped alignments
		if strings.Index(rec.Cigar.String(), "H") != -1 {
			continue
		}
		if indices, ok := signalingReads[rec.Name]; ok {
			for i := range indices {
				if indices[i].pairNumber == getPairNumber(rec) {
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
