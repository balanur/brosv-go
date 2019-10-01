package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"strconv"
	"strings"

	"github.com/balanur/vcfgo"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

var (
	mode = flag.Int("mode", 3, "Running mode.\n"+
		"1: Generate Signaling Reads\n"+
		"2: Voting\n")
	vcfFile = flag.String("vcf", "", "vcf input file")
	bamFile = flag.String("bam", "", "bam input file")
	refFile = flag.String("ref", "", "reference file")
	sr      = flag.String("sr", "signalingReads40.txt", "Signaling Reads file")
	workdir = flag.String("workdir", "", "Working directory")
	threads = flag.Int("threads", 0, "number of threads to use (0 = auto)")
	margin  = flag.Int("margin", 0, "number of error bp allowed (0 = auto)")
	help    = flag.Bool("help", false, "display help")
)

var svTag, lbpTag, rbpTag, copyTag sam.Tag
var segmentSize, variance int
var leftCIs, rightCIs, copyCIs map[string]int
var ciStore CIStore
var svStore SVStore

func readVcf(fileName string) (SVStore, CIStore) {
	return readVcfFiltered(fileName, "", "")
}

func readVcfFiltered(fileName string, filter string, samplefilter string) (SVStore, CIStore) {
	svStore := NewSVStore()
	ciStore := NewCIStore()

	f, _ := os.Open(fileName)
	rdr, err := vcfgo.NewReader(f, false)
	if err != nil {
		panic(err)
	}

	var filter2 string
	if filter == "tandup" {
		filter = "DUP"
		filter2 = "DUP:TANDEM"
	} else if filter == "intdup" {
		filter = "DUP"
		filter2 = "DUP:ISP"
	}
	SVcount := 0
	for {
		variant := rdr.Read()

		if variant == nil {
			break
		}

		svType, _ := variant.Info().Get("SVTYPE")
		sample, _ := variant.Info().Get("SAMPLE")

		if !strings.Contains(svType.(string), filter) || !strings.Contains(variant.Alt()[0], filter2) {
			continue
		}

		if sample != nil && !strings.Contains(sample.(string), samplefilter) {
			continue
		}

		var tempSV SV
		tempSV.Start = int(variant.Pos)
		endPosition, _ := variant.Info().Get("END")
		tempSV.End = endPosition.(int)
		tempSV.Chromosome = variant.Chromosome

		if tempSV.Chromosome == "MT" {
			continue
		}

		tempSV.Type = svType.(string)
		tempSV.id = strings.TrimSpace(variant.Id())

		svsize := tempSV.End - tempSV.Start

		svStore.add(tempSV)
		SVcount++

		var leftInterval Interval
		if ciPos, err := variant.Info().Get("CIPOS"); err == nil {
			leftInterval.head = tempSV.Start + ciPos.([]int)[0]
			leftInterval.tail = tempSV.Start + ciPos.([]int)[1]
			leftInterval.svId = tempSV.id
			leftInterval.side = leftCI
		} else {
			leftInterval.head = tempSV.Start
			leftInterval.tail = tempSV.Start
			leftInterval.svId = tempSV.id
			leftInterval.side = leftCI
		}

		var rightInterval Interval
		if ciEnd, err := variant.Info().Get("CIEND"); err == nil {
			rightInterval.head = tempSV.End + ciEnd.([]int)[0]
			rightInterval.tail = tempSV.End + ciEnd.([]int)[1]
			rightInterval.svId = tempSV.id
			rightInterval.side = rightCI
		} else {
			rightInterval.head = tempSV.End
			rightInterval.tail = tempSV.End
			rightInterval.svId = tempSV.id
			rightInterval.side = rightCI
		}

		if svsize < segmentSize+100 {
			leftInterval.tail += svsize / 2
			rightInterval.head -= svsize / 2
		} else {
			leftInterval.tail += (segmentSize + 100)
			rightInterval.head -= (segmentSize + 100)
		}
		leftInterval.head -= 100
		rightInterval.tail += 100

		if leftInterval.head < 0 {
			leftInterval.head = 1
		}
		if rightInterval.head < 0 {
			rightInterval.head = 1
		}
		ciStore.add(svStore, leftInterval)
		ciStore.add(svStore, rightInterval)

		var copyInterval Interval
		if filter2 == "DUP:ISP" {
			if strings.Contains(variant.Alt()[0], "DUP:ISP") {
				start, _ := variant.Info().Get("POS2")
				tempSV.copyPos, _ = strconv.Atoi(start.(string))
				copyInterval.head = tempSV.copyPos - (segmentSize + 100)
				copyInterval.tail = tempSV.copyPos + (segmentSize + 100)
				copyInterval.svId = tempSV.id
				copyInterval.side = copyCI

				if copyInterval.head < 0 {
					copyInterval.head = 1
				}

				ciStore.add(svStore, copyInterval)
			}
		}

	}
	fmt.Printf("Number of CIs / SVs %d / %d\n", len(ciStore.ciList), len(svStore.svMap))
	fmt.Printf("Total SV count: %d\n", SVcount)
	return svStore, ciStore
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
	f.Close()
	f, _ = os.Open(bamFilePath)
	defer f.Close()

	bamReader, _ = bam.NewReader(f, *threads)
	defer bamReader.Close()

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

func linkLeftRightCIs(svStore SVStore, ciStore CIStore) {
	leftCIs = make(map[string]int)
	rightCIs = make(map[string]int)
	copyCIs = make(map[string]int)

	for i, interval := range ciStore.ciList {
		if interval.side == leftCI {
			leftCIs[interval.svId] = i
		} else if interval.side == rightCI {
			rightCIs[interval.svId] = i
		} else {
			copyCIs[interval.svId] = i
		}
	}
}

// Organizer functions for each step of the workflow
func extractSignalingReadsMode(svType SVType) {
	fmt.Printf("Running in mode 1 - Signaling read extraction\n")
	extractSignalingInCI(*bamFile, path.Join(*workdir, "cluster.bam"), ciStore, svType)
	setBreakpointTags(path.Join(*workdir, "cluster.bam"), *workdir, ciStore)
}

func votingMode(strType string, sample string) {
	fmt.Printf("Running in mode 2 - Breakpoint Voting \n")
	cmd := exec.Command("samtools", "sort", "-t", "SV", path.Join(*workdir, "cluster_withbp.bam"), "-o", path.Join(*workdir, "sorted.bam"))
	cmd.Run()
	calculateSplitReadSupport(path.Join(*workdir, "sorted.bam"), path.Join(*workdir, "votes.txt"), ciStore, svStore)
	writeRefinedVcf(path.Join(*workdir, "votes.txt"), path.Join(*workdir, "refined.vcf"), *refFile, ciStore, svStore)
	compareWithTruth(path.Join(*workdir, "refined.vcf"), "data/simu/del_true_all.bed", strType, sample, ciStore)

}

func main() {
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}

	svTag = sam.NewTag("SV")
	lbpTag = sam.NewTag("LBP")
	rbpTag = sam.NewTag("RBP")
	copyTag = sam.NewTag("CPY")

	svType := del

	var strType string
	if svType == del {
		strType = "DEL"
	} else if svType == inv {
		strType = "INV"
	} else if svType == ins {
		strType = "INS"
	} else if svType == tandup {
		strType = "tandup"
	} else if svType == intdup {
		strType = "intdup"
	}

	segmentSize, variance = findAverageSegmentSize(*bamFile, 1000, 1000000)
	fmt.Printf("Segment size = %d  Variance = %d\n", segmentSize, variance)
	if svType == all {
		svStore, ciStore = readVcf(*vcfFile)
	} else {
		svStore, ciStore = readVcfFiltered(*vcfFile, strType, "")
	}
	linkLeftRightCIs(svStore, ciStore)

	/*
		writeCIstobed(path.Join(*workdir, "cifile.csv"), ciStore, strType)
		return
		compareWithTruth("data/simu/lumpy_30x.vcf", "data/simu/del_true_all.bed", "del", "", ciStore)
		return
	*/

	switch *mode {
	case 1:
		extractSignalingReadsMode(svType)
	case 2:
		votingMode(strType, "")
	case 3:
		alignmentMode()
	default:
		extractSignalingReadsMode(svType)
		votingMode(strType, "")
	}

}

/*
	./brosv-go -vcf data/tardis_40x.vcf -bam data/cnv_1200_40x.bam -ref data/human_g1k_v37.fasta -threads 8 -margin 5 -mode 2 -workdir dels/
*/
