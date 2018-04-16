package main

import (
	"bufio"
	"os"
	"strconv"

	"github.com/biogo/hts/sam"
)

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

func getPairNumber(rec *sam.Record) int {
	if rec.Flags&sam.Read1 != 0 {
		return 1
	}
	return 2
}

func getMappingOri(rec *sam.Record) int {
	if rec.Flags&sam.Reverse == 0 {
		return 1
	}
	return 2
}
