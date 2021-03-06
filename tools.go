package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/sam"
)

func getSecondary(rec *sam.Record) int {
	if rec.Flags&sam.Secondary != 0 {
		return 1
	}
	return 0
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

func writeCIstobed(cifile string, ciStore CIStore, strType string) {

	//Output file
	g, _ := os.Create(cifile)
	defer g.Close()
	writer := bufio.NewWriter(g)

	for id, i := range leftCIs {
		j := rightCIs[id]
		k := copyCIs[id]
		left := ciStore.ciList[i]
		right := ciStore.ciList[j]
		copy := ciStore.ciList[k]
		if left.head >= left.tail {
			fmt.Printf("wtf %d %d\n", left.head, left.tail)
		}
		if right.head >= right.tail {
			fmt.Printf("wtf2 %d %d\n", right.head, right.tail)
		}
		//writer.WriteString(svStore.svMap[id].Chromosome + "\t" + strconv.Itoa(left.head) + "\t" + strconv.Itoa(left.tail) + "\n")
		writer.WriteString(svStore.svMap[id].Chromosome + "," + strconv.Itoa(left.tail-left.head) + "," + strconv.Itoa(right.tail-right.head) + "\n")
		//writer.WriteString(svStore.svMap[id].Chromosome + "\t" + strconv.Itoa(right.head) + "\t" + strconv.Itoa(right.tail) + "\n")
		if strType == "intdup" {
			if copy.head >= copy.tail {
				fmt.Printf("wtf3 %d %d\n", copy.head, copy.tail)
			}
			writer.WriteString(svStore.svMap[id].Chromosome + "\t" + strconv.Itoa(copy.head) + "\t" + strconv.Itoa(copy.tail) + "\n")
		}
	}
	writer.Flush()
}

func simStatistics(simfile string) {
	f, _ := os.Open(simfile)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	indel := 0
	sv1 := 0
	sv2 := 0
	sv3 := 0
	sv4 := 0

	indel_i := 0
	sv1_i := 0
	sv2_i := 0
	sv3_i := 0
	sv4_i := 0

	inter := 0
	invert := 0

	if strings.Contains(simfile, "bed") {

		for scanner.Scan() {
			words := strings.Fields(scanner.Text())
			e, _ := strconv.Atoi(words[2])
			s, _ := strconv.Atoi(words[1])
			len := e - s
			if len <= 50 {
				indel++
			} else if len <= 500 {
				sv1++
			} else if len <= 5000 {
				sv2++
			} else if len <= 10000 {
				sv3++
			} else {
				sv4++
			}
		}
	} else {

		for scanner.Scan() {
			words := strings.Fields(scanner.Text())
			len, _ := strconv.Atoi(words[3])
			duptype := words[11]
			if duptype == "tandem" {
				if len <= 50 {
					indel++
				} else if len <= 500 {
					sv1++
				} else if len <= 5000 {
					sv2++
				} else if len <= 10000 {
					sv3++
				} else {
					sv4++
				}
			} else if duptype == "interspersed" || duptype == "inverted" {
				if len <= 50 {
					indel_i++
				} else if len <= 500 {
					sv1_i++
				} else if len <= 5000 {
					sv2_i++
				} else if len <= 10000 {
					sv3_i++
				} else {
					sv4_i++
				}

				if duptype == "interspersed" {
					inter++
				}
				if duptype == "inverted" {
					invert++
				}
			}
		}
	}
	fmt.Printf("Tandem\n")
	fmt.Printf("1-50:\t%d\n50-500:\t%d\n500-5000:\t%d\n5000-10000:\t%d\n>10000:\t%d\n", indel, sv1, sv2, sv3, sv4)
	fmt.Printf("Interspersed\n")
	fmt.Printf("1-50:\t%d\n50-500:\t%d\n500-5000:\t%d\n5000-10000:\t%d\n>10000:\t%d\n", indel_i, sv1_i, sv2_i, sv3_i, sv4_i)
	fmt.Printf("interspersed: %d\ninverted: %d\n", invert, invert)
}
