package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strconv"

	"github.com/biogo/hts/bam"
)

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
