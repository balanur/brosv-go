package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
)

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
			loc = auxValue(rec.AuxFields.Get(rbpTag))
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

func combineResults(voteFile string /*contigFile string,*/, outfilePath string, ciStore CIStore, svStore SVStore) {
	f, _ := os.Open(voteFile)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	/*f2, _ := os.Open(contigFile)
	defer f2.Close()
	bamReader, _ := bam.NewReader(f2, *threads)
	defer bamReader.Close()
	*/
	f3, _ := os.Create(outfilePath)
	defer f3.Close()

	writer := bufio.NewWriter(f3)

	f4, _ := os.Open("data/tardis_40x.vcf")
	defer f4.Close()
	vcfscanner := bufio.NewScanner(f4)

	leftbp := make(map[string]Loc)
	rightbp := make(map[string]Loc)
	//cigarsup := make(map[string][]int)

	// split read support
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		if words[0][0] == 'c' {
			ciId, _ := strconv.Atoi(words[1])
			svId := ciStore.ciList[ciId].svId
			side, _ := strconv.Atoi(words[2])
			scanner.Scan()
			words := strings.Fields(scanner.Text())
			pos, _ := strconv.Atoi(words[0])
			support, _ := strconv.Atoi(words[1])
			// fill sv maps
			if side == 1 {
				leftbp[svId] = Loc{Pos: pos, VoteNum: support}
			} else {
				rightbp[svId] = Loc{Pos: pos, VoteNum: support}
			}
		}
	} /*
		// contig support
		for {
			rec, err := bamReader.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Fatalf("error reading bam: %v", err)
			}
			if isValidSplit(rec.Cigar) {
				ciId := auxValue(rec.AuxFields.Get(svTag))
				lbp := auxValue(rec.AuxFields.Get(lbpTag))
				rbp := auxValue(rec.AuxFields.Get(rbpTag))
				svId := ciStore.ciList[ciId].svId
				cigarsup[svId] = append(cigarsup[svId], lbp)
				cigarsup[svId] = append(cigarsup[svId], rbp)
			}
		}*/

	// write map to vcf
	for vcfscanner.Scan() {
		if vcfscanner.Text()[0:2] == "##" {
			writer.WriteString(vcfscanner.Text() + "\n")
		}
	}
	// fix this wtf is this
	writer.WriteString("##INFO=<ID=SRSUPL,Number=1,Type=Integer,Description=\"Number of supporting split reads on left\">\n")
	writer.WriteString("##INFO=<ID=SRSUPR,Number=1,Type=Integer,Description=\"Number of supporting split reads on right\">\n")
	writer.WriteString("##INFO=<ID=CONTIGSUPLBP,Number=1,Type=Integer,Description=\"Pos of left breakpoint if supported by assembly\">\n")
	writer.WriteString("##INFO=<ID=CONTIGSUPRBP,Number=1,Type=Integer,Description=\"Pos of right breakpoint if supported by assembly\">\n")
	writer.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcnv_1000_ref\n")
	for svId := range leftbp {
		_sv := svStore.get(svId)
		writer.WriteString(_sv.Chromosome + "\t" + strconv.Itoa(leftbp[svId].Pos) + "\t" + _sv.id + "\t.\t" + _sv.Type + "\t255\tPASS\t")
		svlen := rightbp[svId].Pos - leftbp[svId].Pos
		writer.WriteString("END=" + strconv.Itoa(rightbp[svId].Pos) + ";SVLEN=" + strconv.Itoa(svlen))
		writer.WriteString(";SRSUPL=" + strconv.Itoa(leftbp[svId].VoteNum) + ";SRSUPR=" + strconv.Itoa(rightbp[svId].VoteNum))
		//if len(cigarsup[svId]) != 0 {
		//	writer.WriteString(";CONTIGSUPLBP=" + strconv.Itoa(cigarsup[svId][0]) + ";CONTIGSUPRBP=" + strconv.Itoa(cigarsup[svId][1]) + "\n")
		//} else {
		writer.WriteString("\n")
		//}
	}
	writer.Flush()
}

func compareWithTruth(resultfile string, truthfile string) {
	f, _ := os.Open(resultfile)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	f2, _ := os.Open(truthfile)
	defer f2.Close()
	scanner2 := bufio.NewScanner(f2)

	var truth []SV
	var result []SV

	count := 0
	for scanner2.Scan() {
		words := strings.Fields(scanner2.Text())
		_start, _ := strconv.Atoi(words[1])
		_end, _ := strconv.Atoi(words[2])
		_len := _end - _start //strconv.Atoi(words[3])
		_type := words[5]
		if _type == "del" && _len <= 10000 {
			truth = append(truth, SV{id: ".", Chromosome: words[0][3:], Start: _start, End: _end, Type: _type})
			count++
		}
	}
	sort.Slice(truth, func(i, j int) bool { return truth[i].Start < truth[j].Start })

	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		if words[0][0] == '#' {
			continue
		}
		_start, _ := strconv.Atoi(words[1])
		info := words[7]
		pos := strings.Index(info, ";")
		_end, _ := strconv.Atoi(info[4:pos])
		result = append(result, SV{id: words[2], Chromosome: words[0], Start: _start, End: _end, Type: "del"})
	}
	sort.Slice(result, func(i, j int) bool { return result[i].Start < result[j].Start })

	fmt.Printf("result len %d \n", len(result))
	// compare
	margin := 5
	j := 0
	lTRUE := 0
	rTRUE := 0
	for i := 0; i < len(result) && j < len(truth); {

		if math.Abs(float64(result[i].Start-truth[j].Start)) == float64(margin) {
			fmt.Printf("%s L %d\n", result[i].id, result[i].Start)
		}
		if math.Abs(float64(result[i].End-truth[j].End)) == float64(margin) {
			fmt.Printf("%s R %d\n", result[i].id, result[i].End)
		}

		if result[i].Start-truth[j].Start > margin {
			j++
		} else if truth[j].Start-result[i].Start > margin {
			i++
		} else {
			lTRUE++
			i++
			j++
		}
	}
	j = 0
	for i := 0; i < len(result) && j < len(truth); {
		if result[i].End-truth[j].End > margin {
			j++
		} else if truth[j].End-result[i].End > margin {
			i++
		} else {
			rTRUE++
			i++
			j++
		}
	}
	fmt.Printf("LEFT bp found %d\nRIGHT bp found %d\n", lTRUE, rTRUE)
}
