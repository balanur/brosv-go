package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
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
		ref.SetID(id)
		seenRefs[words[0]] = ref
	}

	return seenRefs
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
func getMateOri(rec *sam.Record) int {
	if rec.Flags&sam.MateReverse == 0 {
		return 1
	}
	return 2
}

// ParseFaiLine go
func ParseFaiLine(line string) FaiEntry {
	lineScanner := bufio.NewScanner(strings.NewReader(line))
	lineScanner.Split(bufio.ScanWords)
	var counter int
	var result FaiEntry
	for lineScanner.Scan() {
		word := lineScanner.Text()
		var e error
		switch counter {
		case 0:
			result.title = word
		case 1:
			result.length, e = strconv.ParseInt(word, 10, 64)
		case 2:
			result.offset, e = strconv.ParseInt(word, 10, 64)
		}
		if e != nil {
			log.Fatal(e)
		}
		counter++
	}
	return result
}

// ReadChr go
func ReadChr(file *os.File, entry FaiEntry) Chromosome {
	var result Chromosome
	result.title = entry.title

	var buffer bytes.Buffer
	fmt.Printf("Seeking to %d for %s\n", entry.offset, entry.title)
	file.Seek(entry.offset, 0)
	scanner := bufio.NewScanner(file)
	scanner.Split(bufio.ScanLines)

	for scanner.Scan() {
		line := scanner.Text()
		if len(line) > 0 && line[0] == '>' {
			result.content = buffer.String()
			buffer.Reset()
			break
		} else {
			buffer.WriteString(line)
		}
	}
	return result
}

func readReference(referencePath string) Genome {
	var genome Genome
	genome.chmMap = make(map[string]int)
	faiFile, err := os.Open(referencePath + ".fai")
	if err != nil {
		log.Fatal(err)
	}
	defer faiFile.Close()

	fastaFile, err := os.Open(referencePath)
	if err != nil {
		log.Fatal(err)
	}
	defer fastaFile.Close()

	log.Println("Reading reference genome from", referencePath)
	scanner := bufio.NewScanner(faiFile)
	for scanner.Scan() {
		entry := ParseFaiLine(scanner.Text())
		genome.faiEntries = append(genome.faiEntries, entry)

		chm := ReadChr(fastaFile, entry)
		genome.chms = append(genome.chms, chm)
		genome.chmMap[chm.title] = len(genome.chms) - 1
		fmt.Printf("Loaded %s\t\t\t\r", chm.title)
	}
	return genome
}

func auxValue(aux sam.Aux) int {
	if aux == nil {
		fmt.Printf("aux is nil\n")
	}
	switch aux.Value().(type) {
	case int8:
		return int(aux.Value().(int8))
	case int16:
		return int(aux.Value().(int16))
	case int32:
		return int(aux.Value().(int32))
	case int64:
		return int(aux.Value().(int64))
	case uint8:
		return int(aux.Value().(uint8))
	case uint16:
		return int(aux.Value().(uint16))
	case uint32:
		return int(aux.Value().(uint32))
	case uint64:
		return int(aux.Value().(uint64))
	default:
		return aux.Value().(int)
	}
}

func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func Complement(s string) string {
	runes := []rune(s)
	for i := 0; i < len(runes); i++ {
		switch runes[i] {
		case 'A':
			runes[i] = 'T'
		case 'T':
			runes[i] = 'A'
		case 'G':
			runes[i] = 'C'
		case 'C':
			runes[i] = 'G'
		default:
			runes[i] = 'N'
		}
	}
	return string(runes)
}

func computeCIGAR(aL string, bL string, cL string) (sam.Cigar, int) {
	M := 0
	D := 0
	I := 0
	S := 0
	var cigar sam.Cigar
	ciglen := 0

	for i := 0; i < len(aL); i++ {
		// Soft clip
		if cL[i] == 'S' {
			S++
			if I != 0 {
				ciglen += I
				cigar = append(cigar, sam.NewCigarOp(sam.CigarInsertion, I))
				I = 0
			}
			if D != 0 {
				cigar = append(cigar, sam.NewCigarOp(sam.CigarDeletion, D))
				D = 0
			}
			if M != 0 {
				ciglen += M
				cigar = append(cigar, sam.NewCigarOp(sam.CigarMatch, M))
				M = 0
			}
		} else if cL[i] == '|' { // Alignment Match
			M++
			if I != 0 {
				ciglen += I
				cigar = append(cigar, sam.NewCigarOp(sam.CigarInsertion, I))
				I = 0
			}
			if S != 0 {
				ciglen += S
				cigar = append(cigar, sam.NewCigarOp(sam.CigarSoftClipped, S))
				S = 0
			}
			if D != 0 {
				cigar = append(cigar, sam.NewCigarOp(sam.CigarDeletion, D))
				D = 0
			}
		} else if cL[i] == ' ' && aL[i] == '-' { // Insention to reference
			I++
			if M != 0 {
				ciglen += M
				cigar = append(cigar, sam.NewCigarOp(sam.CigarMatch, M))
				M = 0
			}
			if S != 0 {
				ciglen += S
				cigar = append(cigar, sam.NewCigarOp(sam.CigarSoftClipped, S))
				S = 0
			}
			if D != 0 {
				cigar = append(cigar, sam.NewCigarOp(sam.CigarDeletion, D))
				D = 0
			}
		} else if cL[i] == ' ' && bL[i] == '-' { // Deletion from reference
			D++
			if I != 0 {
				ciglen += I
				cigar = append(cigar, sam.NewCigarOp(sam.CigarInsertion, I))
				I = 0
			}
			if M != 0 {
				ciglen += M
				cigar = append(cigar, sam.NewCigarOp(sam.CigarMatch, M))
				M = 0
			}
			if S != 0 {
				ciglen += S
				cigar = append(cigar, sam.NewCigarOp(sam.CigarSoftClipped, S))
				S = 0
			}
		} else { // Mismatch
			M++
			if I != 0 {
				ciglen += I
				cigar = append(cigar, sam.NewCigarOp(sam.CigarInsertion, I))
				I = 0
			}
			if S != 0 {
				ciglen += S
				cigar = append(cigar, sam.NewCigarOp(sam.CigarSoftClipped, S))
				S = 0
			}
			if D != 0 {
				cigar = append(cigar, sam.NewCigarOp(sam.CigarDeletion, D))
				D = 0
			}
		}
	}
	// last sequence
	if I != 0 {
		ciglen += I
		cigar = append(cigar, sam.NewCigarOp(sam.CigarInsertion, I))
		I = 0
	}
	if M != 0 {
		ciglen += M
		cigar = append(cigar, sam.NewCigarOp(sam.CigarMatch, M))
		M = 0
	}
	if S != 0 {
		ciglen += S
		cigar = append(cigar, sam.NewCigarOp(sam.CigarSoftClipped, S))
		S = 0
	}
	if D != 0 {
		cigar = append(cigar, sam.NewCigarOp(sam.CigarDeletion, D))
		D = 0
	}
	return cigar, ciglen
}

func checkSplits(bamFilePath string) {
	fmt.Printf("Checking bam\n")
	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	cnt := 0
	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}

		read := string(rec.Seq.Expand())
		if len(read) < 60 {
			// reads with hard clips - dont know what to do with them
			fmt.Printf("small read %s with %d %d %s\n", rec.Cigar.String(), rec.RefID(), rec.Pos, read)
		}
		cnt++
	}
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

// brand new function to vote breakpoint locations
func calculateSplitReadSupport(bamFilePath string, outfile string, ciStore CIStore, svStore SVStore) {

	//Input file
	f, _ := os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %v:\n", err)
	}

	//Output file
	g, _ := os.Create(outfile)
	defer g.Close()
	writer := bufio.NewWriter(g)

	breakpoints := make(map[int]map[int]int)
	current := -1
	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		// get ci index of read
		ciIndex := auxValue(rec.AuxFields.Get(svTag))
		// update if you pass to a new ci
		if current != ciIndex {
			current = ciIndex
			breakpoints[ciIndex] = make(map[int]int)
		}

		// get bp loc left or right
		var loc int
		if ciStore.ciList[ciIndex].side == leftCI {
			loc = auxValue(rec.AuxFields.Get(lbpTag))
		} else {
			loc = auxValue(rec.AuxFields.Get(rbpTag)) - 1
		}

		// check if split provides enough support (not like 97M-3M, 100M etc.)
		// and len > 80
		if isValidSplit(rec.Cigar) {
			// update num of votes
			if _, exist := breakpoints[current][loc]; exist {
				breakpoints[current][loc]++
			} else {
				breakpoints[current][loc] = 1
			}
		}
	}

	// write result to file
	for k, v := range breakpoints {
		var side int
		if ciStore.ciList[k].side == leftCI {
			side = 1
		} else {
			side = 2
		}

		var list []Loc
		for pos, votes := range v {
			list = append(list, Loc{Pos: pos, VoteNum: votes})
		}

		// if there is no support dont write it
		if len(list) > 0 {
			writer.WriteString("ci " + strconv.Itoa(k) + " " + strconv.Itoa(side) + "\n")
		}

		sort.Slice(list, func(i, j int) bool { return list[i].VoteNum > list[j].VoteNum })
		for _, val := range list {
			writer.WriteString(strconv.Itoa(val.Pos) + " " + strconv.Itoa(val.VoteNum) + "\n")
		}
	}
	writer.Flush()
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
		//currentSV := svStore.get(currentCI.svId)
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
func formatQual(q []byte) []byte {
	for _, v := range q {
		if v != 0xff {
			a := make([]byte, len(q))
			for i, p := range q {
				a[i] = p + 33
			}
			return a
		}
	}
	return []byte{'*'}
}
