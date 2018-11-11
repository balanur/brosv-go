package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/brentp/vcfgo"
)

var (
	mode = flag.Int("mode", 4, "Running mode.\n"+
		"0: Generate Signaling Reads\n"+
		"1: Generate Clusters\n"+
		"2: Assembly\n"+
		"3: Alignment\n")
	vcfFile = flag.String("vcf", "", "vcf input file")
	bamFile = flag.String("bam", "", "bam input file")
	refFile = flag.String("ref", "", "reference file")
	sr      = flag.String("sr", "signalingReads40.txt", "Signaling Reads file")
	workdir = flag.String("workdir", "", "Working directory")
	threads = flag.Int("threads", 0, "number of threads to use (0 = auto)")
	help    = flag.Bool("help", false, "display help")
)

var svTag, lbpTag, rbpTag sam.Tag
var segmentSize, variance int
var leftCIs, rightCIs map[string]int
var ciStore CIStore
var svStore SVStore

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

	SVcount := 0
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

		delsize := tempSV.End - tempSV.Start
		if delsize > 10000 {
			continue
		}
		svStore.add(tempSV)
		SVcount++

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

		// All SVs (exact or not) taken
		leftInterval.head -= (segmentSize + 100)
		ciStore.add(svStore, leftInterval)
		rightInterval.tail += (segmentSize + 100)
		ciStore.add(svStore, rightInterval)

	}
	fmt.Printf("Number of CIs / SVs %d / %d\n", len(ciStore.ciList), len(svStore.svMap))
	fmt.Printf("Total dels count: %d\n", SVcount)
	return svStore, ciStore
}

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
				} else if svType == inv {
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

func findAverageSegmentSize(bamFilePath string, thr int, N int) (int, int) {
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
	mean := sum / i
	i = 0
	sum = 0
	for {
		rec, err := bamReader.Read()
		if err == io.EOF || i >= N {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}
		if rec.Flags&sam.Paired != 0 && rec.Flags&sam.ProperPair != 0 && rec.TempLen > 0 && rec.TempLen < thr {
			sum += (rec.TempLen - mean) * (rec.TempLen - mean)
			i++
		}
	}
	variance := int(math.Sqrt(float64(sum / i)))
	return mean, variance
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

		if indices, ok := signalingReads[rec.Name]; ok {
			for i := range indices {
				if indices[i].pairNumber == getPairNumber(rec) {
					//interval := ciStore.get(indices[i].ciIndex)

					/*if interval.side == leftCI && indices[i].mappingOri == 2 {
						continue
					}
					if interval.side == rightCI && indices[i].mappingOri == 1 {
						continue
					}*/

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
func alignContigs(refFilePath string, clusterBamPath string, contigfile string, outfilePath string, ciStore CIStore, svStore SVStore, seenRefs map[string]*sam.Reference) {
	f, _ := os.Open(contigfile)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	f2, _ := os.Open(clusterBamPath)
	defer f.Close()
	bamReader, _ := bam.NewReader(f2, *threads)
	defer bamReader.Close()

	g, _ := os.Create(outfilePath + ".bam")
	//g, _ := os.Create("temp/test.bam")
	defer g.Close()
	alWriter, _ := bam.NewWriter(g, bamReader.Header(), 0)
	defer alWriter.Close()

	ref := readReference(refFilePath)
	var refL, refR string
	var currentCI, l, r Interval
	var currentSV SV
	var ciIndex int
	count := 0
	allcount := 0

	for scanner.Scan() {
		line := scanner.Text()
		// processing new CI
		if line[0] == 'c' {
			words := strings.Fields(scanner.Text())
			ciIndex, _ = strconv.Atoi(words[1])
			fmt.Printf("Processing ci %d reads mapped so far %d / %d\n", ciIndex, count, allcount)
			currentCI = ciStore.ciList[ciIndex]
			//currentSV = svStore.get(currentCI.svId)
			if currentCI.side == leftCI {
				l = currentCI
				r = ciStore.ciList[rightCIs[currentCI.svId]]
			} else {
				l = ciStore.ciList[leftCIs[currentCI.svId]]
				r = currentCI
			}
			refL = ref.getChr(currentSV.Chromosome).content[l.head : l.tail+1]
			refR = ref.getChr(currentSV.Chromosome).content[r.head : r.tail+1]
		} else if line[0] == '>' {
			continue
		} else {
			allcount++
			contig := scanner.Text()
			var result AlignmentResult
			var flags sam.Flags
			initAligner(len(contig))
			// Reverse compelement the left side
			if currentCI.side == leftCI {
				contig = Complement(Reverse(contig))
				result = align(refL, refR, contig)
				flags = sam.Reverse
			} else {
				result = align(refL, refR, contig)
				flags = sam.MateReverse
			}

			if (result.identityL >= 0.95 && result.identityR >= 0.95) || (result.identityL == -1 && result.identityR >= 0.95) || (result.identityL >= 0.95 && result.identityR == -1) {
				cigarL, _ := computeCIGAR(result.aL, result.bL, result.cL)
				cigarR, _ := computeCIGAR(result.aR, result.bR, result.cR)

				var cigar sam.Cigar
				dellen := (r.head + result.rbp) - (l.head + result.lbp) - 1
				if dellen > 0 {
					delreg := sam.NewCigarOp(sam.CigarDeletion, dellen)
					cigar = append(cigarL, delreg)
					cigar = append(cigar, cigarR...)
				} else {
					cigar = append(cigarL, cigarR...)
				}
				chr := ref.getChr(currentSV.Chromosome)
				ref, _ := sam.NewReference(chr.title, "", "", len(chr.content), nil, nil)
				ref.SetID(seenRefs[ref.Name()].ID())

				qual := make([]byte, len(contig))
				for i := range qual {
					qual[i] = 34
				}
				rec, err := sam.NewRecord(strconv.Itoa(count), ref, nil, result.pos+l.head, -1, len(contig), 0, cigar, []byte(contig), qual, []sam.Aux{})
				if err != nil {
					log.Fatalf("problem %v\n", err)
				}
				rec.Flags = flags
				// lbp & rbp tags
				newAux, _ := sam.NewAux(lbpTag, result.lbp+l.head)
				rec.AuxFields = append(rec.AuxFields, newAux)
				newAux2, _ := sam.NewAux(rbpTag, result.rbp+r.head)
				rec.AuxFields = append(rec.AuxFields, newAux2)
				// sv tag
				newAux3, _ := sam.NewAux(svTag, ciIndex)
				rec.AuxFields = append(rec.AuxFields, newAux3)
				err = alWriter.Write(rec)
				count++
			}
		}
	}
	fmt.Printf("All contigs vs. mapped ones: %d / %d\n", allcount, count)
}

func alignClusters(refFilePath string, clusterBamPath string, outfilePath string, svStore SVStore, ciStore CIStore) {
	f, _ := os.Open(clusterBamPath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	if ok, err := bgzf.HasEOF(f); err != nil || !ok {
		log.Fatalf("could not open file %q:", err)
	}

	f2, _ := os.Create(outfilePath + ".txt")
	defer f2.Close()

	writer := bufio.NewWriter(f2)

	g1, _ := os.Create(outfilePath + ".bam")
	defer g1.Close()

	alWriter, _ := bam.NewWriter(g1, bamReader.Header(), 0)
	defer alWriter.Close()

	ref := readReference(refFilePath)
	var refL, refR string
	var currentCI, l, r Interval
	var currentSV SV
	ciIndex := -1
	var side, before int
	count := 0
	allcount := 0

	for {
		rec, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}
		temp := auxValue(rec.AuxFields.Get(svTag))
		if ciIndex != temp {
			ciIndex = temp
			fmt.Printf("Processing ci %d\n", ciIndex)
			currentCI = ciStore.ciList[ciIndex]
			currentSV = svStore.get(currentCI.svId)
			if currentCI.side == leftCI {
				side = 1
				l = currentCI
				r = ciStore.ciList[rightCIs[currentCI.svId]]
			} else {
				side = 2
				l = ciStore.ciList[leftCIs[currentCI.svId]]
				r = currentCI
			}
			refL = ref.getChr(currentSV.Chromosome).content[l.head : l.tail+1]
			refR = ref.getChr(currentSV.Chromosome).content[r.head : r.tail+1]
		}
		read := string(rec.Seq.Expand())
		/*if rec.Flags&sam.Reverse == 0 && currentCI.side == leftCI {
			read = Complement(Reverse(read))
		}
		if rec.Flags&sam.Reverse != 0 && currentCI.side == rightCI {
			read = Complement(Reverse(read))
		}*/

		var result AlignmentResult
		initAligner(len(read))
		if currentCI.side == leftCI {
			result = align(refL, refR, read)
			rec.Flags = sam.Paired + sam.ProperPair + sam.Reverse
		} else {
			result = align(refL, refR, read)
			rec.Flags = sam.Paired + sam.ProperPair + sam.MateReverse
		}

		if (result.identityL >= 0.95 && result.identityR >= 0.95) || (result.identityL == -1 && result.identityR >= 0.95) || (result.identityL >= 0.95 && result.identityR == -1) {
			if before != ciIndex {
				writer.WriteString(strconv.Itoa(ciIndex) + " " + strconv.Itoa(side) + ":\n")
				before = ciIndex
			}
			writer.WriteString("chr:" + currentSV.Chromosome + " lbp: " + strconv.Itoa(result.lbp+l.head) + " rbp: " + strconv.Itoa(result.rbp+r.head) + "\n")
			writer.WriteString(result.aL + "\n" + result.cL + "\n" + result.bL + "\nidentityL: " + strconv.FormatFloat(result.identityL, 'f', 3, 64) + "\n\n")
			writer.WriteString(result.aR + "\n" + result.cR + "\n" + result.bR + "\nidentityR: " + strconv.FormatFloat(result.identityR, 'f', 3, 64) + "\n\n")
			cigarL, _ := computeCIGAR(result.aL, result.bL, result.cL)
			cigarR, _ := computeCIGAR(result.aR, result.bR, result.cR)

			var cigar sam.Cigar
			dellen := (r.head + result.rbp) - (l.head + result.lbp) - 1
			if dellen > 0 {
				delreg := sam.NewCigarOp(sam.CigarDeletion, dellen)
				cigar = append(cigarL, delreg)
				cigar = append(cigar, cigarR...)
			} else {
				cigar = append(cigarL, cigarR...)
			}
			rec.Pos = result.pos + l.head
			rec.Ref.SetName(currentSV.Chromosome)
			rec.Cigar = cigar
			newAux, _ := sam.NewAux(lbpTag, result.lbp+l.head)
			rec.AuxFields = append(rec.AuxFields, newAux)
			newAux, _ = sam.NewAux(rbpTag, result.rbp+r.head)
			rec.AuxFields = append(rec.AuxFields, newAux)
			alWriter.Write(rec)
			count++
		}
		allcount++
	}
	writer.Flush()
	fmt.Printf("All signaling reads vs. mapped ones: %d / %d\n", allcount, count)
}

func linkLeftRightCIs(svStore SVStore, ciStore CIStore) {
	leftCIs = make(map[string]int)
	rightCIs = make(map[string]int)

	for i, interval := range ciStore.ciList {
		if interval.side == leftCI {
			leftCIs[interval.svId] = i
		} else {
			rightCIs[interval.svId] = i
		}
	}
}

func assembleReads(clusterBamPath string) {
	f, _ := os.Open(clusterBamPath)
	defer f.Close()

	bamReader, _ := bam.NewReader(f, *threads)
	defer bamReader.Close()

	g, _ := os.Create("tempfile.fastq")
	defer g.Close()
	writer := bufio.NewWriter(g)

	g2, _ := os.Create(*workdir + "simu.contigs")
	defer g2.Close()
	writer2 := bufio.NewWriter(g2)

	num := 0
	left := 0
	right := 0
	currentCI := -1
	for {
		rec, err := bamReader.Read()
		var ci int
		if err != io.EOF {
			ci = auxValue(rec.AuxFields.Get(svTag))
		}

		if currentCI == -1 {
			currentCI = ci
		}
		if currentCI != ci || err == io.EOF {

			if ciStore.ciList[ci].side == leftCI {
				left++
			} else {
				right++
			}
			writer.Flush()
			// do the assembly
			fmt.Printf("Assemling ci %d num of reads %d ", currentCI, num)
			num = 0
			cmd := exec.Command("rm", "-rf", "velv")
			cmd.Run()
			cmd = exec.Command("velveth", "velv/", "60", "-fastq", "-short", "tempfile.fastq", "-strand_specific")
			cmd.Run()
			cmd = exec.Command("velvetg", "velv/", "-unused_reads", "yes", "-min_contig_lgth", "100", "-cov_cutoff", "1")
			cmd.Run()

			os.Remove("tempfile.fastq")
			g, _ = os.Create("tempfile.fastq")
			defer g.Close()
			writer = bufio.NewWriter(g)

			// write results
			contigs, _ := os.Open("velv/contigs.fa")
			defer contigs.Close()
			scanner := bufio.NewScanner(contigs)
			contigcount := 0
			writer2.WriteString("ci " + strconv.Itoa(currentCI))
			for scanner.Scan() {
				line := scanner.Text()
				if line[0] == '>' {
					writer2.WriteString("\n")
					writer2.WriteString(line + "\n")
					contigcount++
				} else {
					writer2.WriteString(line)
				}
			}
			writer2.WriteString("\n")
			fmt.Printf("contig count %d \n", contigcount)

			if err == io.EOF {
				break
			}
			currentCI = ci
		}

		// add read sequence to fasta (read name, seq, ori, qual)
		writer.WriteString("@" + rec.Name + "\n")
		writer.WriteString(string(rec.Seq.Expand()) + "\n")
		writer.WriteString("+\n")
		qual := string(formatQual(rec.Qual))
		writer.WriteString(qual + "\n")
		num++
	}
	writer2.Flush()
	fmt.Printf("left : %d right :  %d\n", left, right)
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
	if mlen+mlen < 50 {
		return false
	}
	return true
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

func extractSignalingReadsMode() {
	fmt.Printf("Running in mode 0 - Signaling read extraction\n")
	var signalingReads map[string][]IndexPair
	signalingReads = matchAlignmentsToSVs(*bamFile, ciStore)
	WriteMapToFile(signalingReads, *workdir+*sr)
}

func clusteringMode() {
	fmt.Printf("Running in mode 1 - Clustering\n")
	fmt.Printf("Clustering for dels < 10,000\n")
	signalingReads := restoreSignalingReads(*workdir + *sr)
	clusters, seenRefs := constructClusterFile(*bamFile, signalingReads, svStore, ciStore, *workdir)
	WriteClusterToFile(clusters, path.Join(*workdir, "clusterSizes.txt"), ciStore, svStore)
	WriteRefsToFile(seenRefs, path.Join(*workdir, "seenRefs.txt"))
	findSVnum(path.Join(*workdir, "clusterSizes.txt"), ciStore)
	cmd := exec.Command("samtools", "sort", "-t", "SV", path.Join(*workdir, "cluster.bam"), "-o", path.Join(*workdir, "cluster_sorted.bam"))
	cmd.Run()
}

func assemblyMode() {
	fmt.Printf("Running in mode 2 - Assembling clusters\n")
	assembleReads(path.Join(*workdir, "cluster_sorted.bam"))
}

func alignmentMode() {
	fmt.Printf("Running in mode 3- Aligning the contigs\n")
	//seenRefs := restoreRefs(path.Join(*workdir, "seenRefs.txt"))
	//alignContigs(*refFile, path.Join(*workdir, "cluster.bam"), path.Join(*workdir, "simu.contigs"), path.Join(*workdir, "alignment_contig"), ciStore, svStore, seenRefs)
	alignClusters(*refFile, *workdir+"cluster_sorted.bam", *workdir+"alignment40", svStore, ciStore)
	extractBreakpointResults(*workdir+"alignment40.bam", *workdir+"supportedSVs.txt", ciStore, svStore)
}

func votingMode() {
	cmd := exec.Command("samtools", "sort", "-t", "SV", path.Join(*workdir, "alignment40_sorted.bam"), "-o", path.Join(*workdir, "sorted.bam"))
	cmd.Run()
	calculateSplitReadSupport(path.Join(*workdir, "sorted.bam"), path.Join(*workdir, "votes.txt"), ciStore, svStore)
	combineResults(path.Join(*workdir, "votes.txt") /*path.Join(*workdir, "alignment_contig.bam"),*/, path.Join(*workdir, "refined.vcf"), ciStore, svStore)
	compareWithTruth(path.Join(*workdir, "refined.vcf"), "data/true_cnv_1200.txt")
}

/*
-vcf /home/asoylev/eccb/tardis_40x.vcf -bam /home/fhormozd/thong/cnvsim/cnv1200-40x/cnv_1200_40x.bam -ref /home/halil/bwa/reference/human_g1k_v37.fasta -threads 8 -mode 3 -workdir extended_delsize/
*/

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	fmt.Printf("This program was infected by Halil. Do not try to resist!\n")

	svTag = sam.NewTag("SV")
	lbpTag = sam.NewTag("LBP")
	rbpTag = sam.NewTag("RBP")

	segmentSize, variance = findAverageSegmentSize(*bamFile, 1000, 1000000)
	fmt.Printf("Segment size = %d  Variance = %d\n", segmentSize, variance)
	svStore, ciStore = readVcfFiltered(*vcfFile, "del")
	linkLeftRightCIs(svStore, ciStore)
	//fmt.Printf("Done reading vcf\nNumber of dels (len < 10,000) %d %d\n", len(leftCIs), len(rightCIs))

	switch *mode {
	case 0:
		// Signaling read extraction
		extractSignalingReadsMode()
	case 1:
		// Clustering & cluster statistics
		clusteringMode()
	case 2:
		// Assembly
		assemblyMode()
	case 3:
		// Alignment
		alignmentMode()
	case 4:
		// Voting and Combining Results
		votingMode()
	default:
		extractSignalingReadsMode()
		clusteringMode()
		assemblyMode()
		alignmentMode()
		votingMode()
	}

}
