package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	//"path"
	"sort"
	"strconv"
	"strings"

	"github.com/balanur/vcfgo"
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
		} else if ciStore.ciList[ciIndex].side == rightCI {
			loc = auxValue(rec.AuxFields.Get(rbpTag))
		} else {
			loc = auxValue(rec.AuxFields.Get(copyTag))
		}

		// update num of votes
		if _, exist := breakpoints[current][loc]; exist {
			breakpoints[current][loc]++
		} else {
			breakpoints[current][loc] = 1
		}

	}

	// write result to file
	for k, v := range breakpoints {
		var side int
		if ciStore.ciList[k].side == leftCI {
			side = 1
		} else if ciStore.ciList[k].side == rightCI {
			side = 2
		} else {
			side = 3
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
		if len(list) < 2 {
			continue
		}
	}
	writer.Flush()

}

func writeRefinedVcf(voteFile string, outfilePath string, refFilePath string, ciStore CIStore, svStore SVStore) {
	f, _ := os.Open(voteFile)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	g, _ := os.Create(outfilePath)
	defer g.Close()

	writer := bufio.NewWriter(g)

	f2, _ := os.Open("data/tardis_40x.vcf")
	defer f2.Close()
	vcfscanner := bufio.NewScanner(f2)

	leftbp := make(map[string]Loc)
	rightbp := make(map[string]Loc)
	copybp := make(map[string]Loc)

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
			} else if side == 2 {
				rightbp[svId] = Loc{Pos: pos, VoteNum: support}
			} else {
				copybp[svId] = Loc{Pos: pos, VoteNum: support}
			}
		}
	}

	// write map to vcf
	for vcfscanner.Scan() {
		if vcfscanner.Text()[0:2] == "##" {
			writer.WriteString(vcfscanner.Text() + "\n")
		}
	}
	// write header
	writer.WriteString("##INFO=<ID=SRSUPL,Number=1,Type=Integer,Description=\"Number of supporting split reads on left\">\n")
	writer.WriteString("##INFO=<ID=SRSUPR,Number=1,Type=Integer,Description=\"Number of supporting split reads on right\">\n")
	writer.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcnv_1000_ref\n")

	ref := readReference(refFilePath)

	for svId := range leftbp {
		_sv := svStore.get(svId)

		// if there is enough support
		if leftbp[svId].VoteNum >= 5 || rightbp[svId].VoteNum >= 5 {

			REF, ALT := getREFALT(ref, _sv, leftbp[svId].Pos-1, rightbp[svId].Pos)
			writer.WriteString(_sv.Chromosome + "\t" + strconv.Itoa(leftbp[svId].Pos) + "\t" + _sv.id + "\t" + REF + "\t" + ALT + "\t255\tPASS\t")
			svlen := rightbp[svId].Pos - leftbp[svId].Pos
			writer.WriteString("SVTYPE=" + _sv.Type + ";END=" + strconv.Itoa(rightbp[svId].Pos) + ";SVLEN=" + strconv.Itoa(svlen))
			if _sv.Type == "DUP:ISP" {
				writer.WriteString(";POS2=" + strconv.Itoa(copybp[svId].Pos))
			}
			writer.WriteString(";SRSUPL=" + strconv.Itoa(leftbp[svId].VoteNum) + ";SRSUPR=" + strconv.Itoa(rightbp[svId].VoteNum))
			if _sv.Type == "DUP:ISP" {
				writer.WriteString(";SRSUPCPY=" + strconv.Itoa(copybp[svId].VoteNum))
			}
			writer.WriteString("\n")
		}
	}
	writer.Flush()
}

func getREFALT(ref Genome, sv SV, start int, end int) (string, string) {
	if end <= start {
		return ".", "."
	}
	if start < 0 {
		return ".", "."
	}

	REF := ref.getChr(sv.Chromosome).content[start:end]

	if sv.Type == "DEL" {
		return REF[0:1], "<DEL>"
	} else if sv.Type == "INV" {
		return REF[0:1], "<INV>"
	} else if sv.Type == "DUP:TANDEM" {
		return REF[0:1], "<DUP:TANDEM>"
	} else if sv.Type == "DUP:ISP" {
		return REF[0:1], "<DUP:ISP>"
	}
	return ".", "."
}

func sortcond(x SV, y SV) bool {
	if x.Chromosome < y.Chromosome {
		return true
	} else if x.Chromosome > y.Chromosome {
		return false
	} else if x.Start < y.Start {
		return true
	} else {
		return false
	}
}

func AbsInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func compareWithTruth(resultfile string, truthfile string, strType string, strSample string, ciStore CIStore) {
	f, _ := os.Open(resultfile)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	var truth []SV
	var result []SV

	var filter string
	if strType == "tandup" {
		strType = "dup"
		filter = "tandem"
	} else if strType == "intdup" {
		strType = "dup"
		filter = "interspersed"
	}

	if strings.Contains(truthfile, "vcf") {
		f2, _ := os.Open(truthfile)
		rdr, err := vcfgo.NewReader(f2, false)
		if err != nil {
			panic(err)
		}
		for {
			variant := rdr.Read()

			if variant == nil {
				break
			}
			svType, _ := variant.Info().Get("SVTYPE")
			_type := svType.(string)
			start := int(variant.Pos)
			endd, _ := variant.Info().Get("END")
			end := endd.(int)
			chr := variant.Chromosome
			truth = append(truth, SV{id: ".", Chromosome: chr, Start: start, End: end, Type: _type})
		}
	} else {
		// get truth file as .bed
		f2, _ := os.Open(truthfile)
		defer f2.Close()
		scanner2 := bufio.NewScanner(f2)

		for scanner2.Scan() {
			words := strings.Fields(scanner2.Text())
			start, _ := strconv.Atoi(words[1])
			end, _ := strconv.Atoi(words[2])
			if len(words) <= 3 {
				truth = append(truth, SV{id: ".", Chromosome: words[0], Start: start, End: end, Type: strType})
				continue
			}
			_type := words[5]

			// not duplication
			if strType == _type && strType != "dup" {
				truth = append(truth, SV{id: ".", Chromosome: words[0][3:], Start: start, End: end, Type: _type})
				continue
			}
			// tandem duplication
			if strType == _type && filter == "tandem" && words[11] == "tandem" {
				truth = append(truth, SV{id: ".", Chromosome: words[0][3:], Start: start, End: end, Type: _type})
				continue
			}
			// interspersed & inverted duplication
			if strType == _type && filter == "interspersed" {
				if words[11] == "interspersed" || words[11] == "inverted" {
					jump, _ := strconv.Atoi(words[10])
					truth = append(truth, SV{id: ".", Chromosome: words[0][3:], Start: start, End: end, Type: _type, copyPos: start + jump})
					continue
				}
			}
		}
	}

	sort.Slice(truth, func(i, j int) bool { return sortcond(truth[i], truth[j]) })
	fmt.Printf("Truth len %d \n", len(truth))

	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		if words[0][0] == '#' {
			continue
		}
		start, _ := strconv.Atoi(words[1])
		info := words[7]
		x := strings.Index(info, "END=")
		y := strings.Index(info[x:len(info)], ";")
		end, _ := strconv.Atoi(info[x+4 : x+y])

		x = strings.Index(info, "SVTYPE=")
		y = strings.Index(info[x:len(info)], ";")
		if y <= 0 {
			y = len(info) - x
		}
		svtype := info[x+7 : x+y]

		if filter == "tandem" {
			strType = "DUP:TANDEM"
		} else if filter == "interspersed" {
			strType = "DUP:ISP"
		}
		if strings.EqualFold(svtype, strType) {

			// interspersed && inverted duplication
			if filter == "interspersed" {
				x := strings.Index(info, "POS2=")
				y := strings.Index(info[x:len(info)], ";")
				copypos, _ := strconv.Atoi(info[x+5 : x+y])
				result = append(result, SV{id: words[2], Chromosome: words[0], Start: start, End: end, Type: strType, copyPos: copypos})
			} else { // all other result
				result = append(result, SV{id: words[2], Chromosome: words[0], Start: start, End: end, Type: strType})
			}
		}
	}
	sort.Slice(result, func(i, j int) bool { return sortcond(result[i], result[j]) })

	fmt.Printf("Result len %d \n", len(result))

	//redundancy check
	var temp2 []SV
	temp2 = append(temp2, result[0])

	for i := 1; i < len(result); i++ {
		if result[i].Chromosome != result[i-1].Chromosome {
			temp2 = append(temp2, result[i])
		} else if result[i-1].End < result[i].Start && result[i-1].Start != result[i].Start && result[i-1].End != result[i].End {
			temp2 = append(temp2, result[i])
		}
	}
	result = temp2
	fmt.Printf("Result len %d (no redundancy)\n", len(result))

	//performanceForCIs(truth, result)
	//return

	var forfig []string
	// compare
	i := 0
	j := 0
	cTRUE := 0
	TPl := 0
	TPr := 0
	FP := 0
	FN := 0
	TP2 := 0
	for i = 0; i < len(result); i++ {
		flag := false
		for j = 0; j < len(truth); j++ {
			if result[i].Chromosome == truth[j].Chromosome {
				if AbsInt(result[i].Start-truth[j].Start) <= *margin {
					forfig = append(forfig, result[i].id)
					TPl++
					flag = true
				}
				if AbsInt(result[i].End-truth[j].End) <= *margin {
					TPr++
					flag = true
				}
			}
		}
		if !flag {
			FP++
		}
	}
	for j = 0; j < len(truth); j++ {
		flag := false
		for i = 0; i < len(result); i++ {
			if result[i].Chromosome == truth[j].Chromosome {
				if AbsInt(result[i].Start-truth[j].Start) <= *margin || AbsInt(result[i].End-truth[j].End) <= *margin {
					TP2++
					flag = true
					break
				}
			}
		}
		if !flag {
			FN++
		}
	}

	fmt.Printf("TPl %d  TPr %d FP %d FN %d\n", TPl, TPr, FP, FN)

	// left breakpoint
	/*for i = 0; i < len(result) && j < len(truth); {
		if result[i].Chromosome > truth[j].Chromosome {
			FNs[truth[j].Start] = truth[j]
			j++
		} else if truth[j].Chromosome > result[i].Chromosome {
			FPs[result[i].id] = result[i]
			i++
		} else {
			if result[i].Start-truth[j].Start > *margin {
				FNs[truth[j].Start] = truth[j]
				j++
			} else if truth[j].Start-result[i].Start > *margin {
				FPs[result[i].id] = result[i]
				i++
			} else {
				verified[result[i].id] = 1
				TPs[result[i].id] = result[i]
				TPs2[truth[j].Start] = truth[j]
				lTRUE++
				i++
				j++
			}
		}
	}

	for ; i < len(result); i++ {
		FPs[result[i].id] = result[i]
	}

	for ; j < len(truth); j++ {
		FNs[truth[j].Start] = truth[j]
	}

	fmt.Printf("TP: %d FP: %d FN: %d\n", len(TPs), len(FPs), len(FNs))

	// right breakpoint
	j = 0
	for i = 0; i < len(result) && j < len(truth); {
		if result[i].Chromosome > truth[j].Chromosome {
			j++
		} else if truth[j].Chromosome > result[i].Chromosome {
			if _, ok := TPs[result[i].id]; !ok {
				FPs[result[i].id] = result[i]
			}
			i++
		} else {
			if result[i].End-truth[j].End > *margin {
				j++
			} else if truth[j].End-result[i].End > *margin {
				if _, ok := TPs[result[i].id]; !ok {
					FPs[result[i].id] = result[i]
				}
				i++
			} else {
				if _, ok := FPs[result[i].id]; ok {
					//fmt.Printf("FP %d %d , %d %d\n", truth[j].Start, truth[j].End, result[i].Start, result[i].End)
					delete(FPs, result[i].id)
				}
				if _, ok := FNs[truth[j].Start]; ok {
					//fmt.Printf("FN %d %d , %d %d\n", truth[j].Start, truth[j].End, result[i].Start, result[i].End)
					delete(FNs, truth[j].Start)
				}
				TPs[result[i].id] = result[i]
				TPs2[truth[j].Start] = truth[j]
				verified[result[i].id] = 1
				rTRUE++
				i++
				j++
			}
		}
	}

	for ; i < len(result); i++ {
		FPs[result[i].id] = result[i]
	}

	for ; j < len(truth); j++ {
		FNs[truth[j].Start] = truth[j]
	}

	fmt.Printf("TP: %d FP: %d FN: %d\n", len(TPs2), len(FPs), len(FNs))
	*/

	if filter == "interspersed" {
		// copy loci for interspersed
		j = 0
		for i = 0; i < len(result) && j < len(truth); {
			if result[i].Chromosome > truth[j].Chromosome {
				j++
			} else if truth[j].Chromosome > result[i].Chromosome {
				i++
			} else {
				if result[i].copyPos-truth[j].copyPos > *margin {
					j++
				} else if truth[j].copyPos-result[i].copyPos > *margin {
					i++
				} else {
					//verified[result[i].id] = 1
					cTRUE++
					i++
				}
			}
		}
	}

	fmt.Printf("Copy bp found %d\n", cTRUE)
}

// tp, fp, fn
func performanceForCIs(truth []SV, result []SV) {

	FPs := make(map[string]SV)
	FNs := make(map[int]SV)
	TPs := make(map[string]SV)
	i := 0
	j := 0

	for i = 0; i < len(result) && j < len(truth); {
		lid := leftCIs[result[i].id]
		if result[i].Chromosome > truth[j].Chromosome {
			FNs[truth[j].Start] = truth[j]
			j++
		} else if truth[j].Chromosome > result[i].Chromosome {
			FPs[result[i].id] = result[i]
			i++
		} else {
			if truth[j].Start < ciStore.ciList[lid].head {
				FNs[truth[j].Start] = truth[j]
				j++
			} else if ciStore.ciList[lid].head <= truth[j].Start && truth[j].Start <= ciStore.ciList[lid].tail {
				TPs[result[i].id] = result[i]
				i++
				j++
			} else if ciStore.ciList[lid].tail < truth[j].Start {
				FPs[result[i].id] = result[i]
				i++
			}
		}
	}
	j = 0
	for i = 0; i < len(result) && j < len(truth); {
		rid := rightCIs[result[i].id]
		if result[i].Chromosome > truth[j].Chromosome {
			_, ok := TPs[result[i].id]
			_, ok2 := FPs[result[i].id]
			if !ok && !ok2 {
				FNs[truth[j].Start] = truth[j]
			}
			j++
		} else if truth[j].Chromosome > result[i].Chromosome {
			if _, ok := TPs[result[i].id]; !ok {
				FPs[result[i].id] = result[i]
			}
			i++
		} else {
			if truth[j].End < ciStore.ciList[rid].head {
				_, ok := TPs[result[i].id]
				_, ok2 := FPs[result[i].id]
				if !ok && !ok2 {
					FNs[truth[j].Start] = truth[j]
				}
				j++
			} else if ciStore.ciList[rid].head <= truth[j].End && truth[j].End <= ciStore.ciList[rid].tail {
				if _, ok := FPs[result[i].id]; ok {
					delete(FPs, result[i].id)
				}
				if _, ok := FNs[truth[j].Start]; ok {
					delete(FNs, truth[j].Start)
				}
				TPs[result[i].id] = result[i]
				i++
				j++
			} else if ciStore.ciList[rid].tail < truth[j].End {
				if _, ok := TPs[result[i].id]; !ok {
					FPs[result[i].id] = result[i]
				}
				i++
			}
		}
	}
	for ; i < len(result); i++ {
		FPs[result[i].id] = result[i]
	}

	for ; j < len(truth); j++ {
		FNs[truth[j].Start] = truth[j]
	}

	fmt.Printf("TP: %d FP: %d FN: %d\n", len(TPs), len(FPs), len(FNs))
}
