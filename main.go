package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	//	"os/exec"
	"path"
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

		svsize := tempSV.End - tempSV.Start
		if svsize > 10000 {
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

	for i, interval := range ciStore.ciList {
		if interval.side == leftCI {
			leftCIs[interval.svId] = i
		} else {
			rightCIs[interval.svId] = i
		}
	}
}

// Organizer functions for each step of the workflow
func extractSignalingReadsMode(svType SVType) {
	fmt.Printf("Running in mode 1 - Signaling read extraction\n")
	//extractSignalingInCI(*bamFile, path.Join(*workdir, "cluster.bam"), ciStore, svType)
	setBreakpointTags(path.Join(*workdir, "cluster.bam"), *workdir, ciStore)
}

func votingMode(strType string) {
	fmt.Printf("Running in mode 2 - Breakpoint Voting \n")
	//cmd := exec.Command("samtools", "sort", "-t", "SV", path.Join(*workdir, "cluster_withbp.bam"), "-o", path.Join(*workdir, "sorted.bam"))
	//cmd.Run()
	//calculateSplitReadSupport(path.Join(*workdir, "sorted.bam"), path.Join(*workdir, "votes.txt"), ciStore, svStore)
	//writeRefinedVcf(path.Join(*workdir, "votes.txt"), path.Join(*workdir, "refined.vcf"), ciStore, svStore)
	compareWithTruth(path.Join(*workdir, "refined.vcf"), "data/true_cnv_1200.txt", strType)
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

	svType := del
	var strType string
	if svType == del {
		strType = "del"
	} else if svType == inv {
		strType = "inv"
	}

	segmentSize, variance = findAverageSegmentSize(*bamFile, 1000, 1000000)
	fmt.Printf("Segment size = %d  Variance = %d\n", segmentSize, variance)
	svStore, ciStore = readVcfFiltered(*vcfFile, strType)
	linkLeftRightCIs(svStore, ciStore)

	switch *mode {
	case 1:
		extractSignalingReadsMode(svType)
	case 2:
		votingMode(strType)
	default:
		extractSignalingReadsMode(svType)
		votingMode(strType)
	}

}

/*
	./brosv-go -vcf data/tardis_40x.vcf -bam data/cnv_1200_40x.bam -ref data/human_g1k_v37.fasta -threads 8 -margin 5 -mode 2 -workdir dels/
*/
