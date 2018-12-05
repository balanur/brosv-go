package main

import (
	"bytes"
	"math"
)

const (
	MATCH_SCORE         = 5
	MISMATCH_SCORE      = -4
	GAP_OPENING_SCORE   = -16
	GAP_EXTENSION_SCORE = -1
)

type AlignmentResult struct {
	aL, cL, bL           string
	aR, bR, cR           string
	identityL, identityR float64
	lbp, rbp, pos        int
}

func align(N int, refL string, refR string, read string, svtype string) AlignmentResult {
	var gapaL, gapbL, gapaR, gapbR, scoreL, scoreR [][]int
	var bestL, bestR []int
	var aL, bL, cL, aR, bR, cR bytes.Buffer

	N++
	MAXSIDE := 100000

	gapaL = make([][]int, MAXSIDE)
	gapbL = make([][]int, MAXSIDE)
	scoreL = make([][]int, MAXSIDE)

	gapaR = make([][]int, MAXSIDE)
	gapbR = make([][]int, MAXSIDE)
	scoreR = make([][]int, MAXSIDE)

	bestL = make([]int, N)
	bestR = make([]int, N)

	for i := 0; i < MAXSIDE; i++ {
		gapaL[i] = make([]int, N)
		gapbL[i] = make([]int, N)
		scoreL[i] = make([]int, N)

		gapaR[i] = make([]int, N)
		gapbR[i] = make([]int, N)
		scoreR[i] = make([]int, N)
	}

	scoreL[0][0], gapaL[0][0], gapbL[0][0], scoreR[0][0], gapaR[0][0], gapbR[0][0] = 0, 0, 0, 0, 0, 0
	bestL[0], bestR[0] = 0, 0

	for i := 1; i < MAXSIDE; i++ {
		val := GAP_OPENING_SCORE + (i-1)*GAP_EXTENSION_SCORE
		scoreL[i][0], scoreR[i][0] = 0, 0
		gapaL[i][0], gapbL[i][0] = val, val
		gapaR[i][0], gapbR[i][0] = val, val

	}
	for i := 1; i < N; i++ {
		val := GAP_OPENING_SCORE + (i-1)*GAP_EXTENSION_SCORE
		scoreL[0][i], scoreR[0][i] = val, val
		gapaL[0][i], gapbL[0][i] = val, val
		gapaR[0][i], gapbR[0][i] = val, val
		bestL[i], bestR[i] = 1, 1
	}

	var result AlignmentResult

	readReversed := Reverse(read)
	if svtype == "del" {
		refR = Reverse(refR)
	}

	for i := 1; i <= len(refL); i++ {
		for j := 1; j <= len(read); j++ {
			gapaL[i][j] = max2(gapaL[i-1][j], scoreL[i-1][j]+GAP_OPENING_SCORE) + GAP_EXTENSION_SCORE
			gapbL[i][j] = max2(gapbL[i][j-1], scoreL[i][j-1]+GAP_OPENING_SCORE) + GAP_EXTENSION_SCORE

			if refL[i-1] == read[j-1] {
				scoreL[i][j] = max3(scoreL[i-1][j-1], gapaL[i-1][j-1], gapbL[i-1][j-1]) + MATCH_SCORE
			} else {
				scoreL[i][j] = max3(scoreL[i-1][j-1], gapaL[i-1][j-1], gapbL[i-1][j-1]) + MISMATCH_SCORE
			}
			// clip of contig until j best mapped to i loc in ref
			if scoreL[i][j] > scoreL[bestL[j]][j] {
				bestL[j] = i
			}
		}
	}

	for i := 1; i <= len(refR); i++ {
		for j := 1; j <= len(read); j++ {
			gapaR[i][j] = max2(gapaR[i-1][j], scoreR[i-1][j]+GAP_OPENING_SCORE) + GAP_EXTENSION_SCORE
			gapbR[i][j] = max2(gapbR[i][j-1], scoreR[i][j-1]+GAP_OPENING_SCORE) + GAP_EXTENSION_SCORE
			if refR[i-1] == readReversed[j-1] {
				scoreR[i][j] = max3(scoreR[i-1][j-1], gapaR[i-1][j-1], gapbR[i-1][j-1]) + MATCH_SCORE
			} else {
				scoreR[i][j] = max3(scoreR[i-1][j-1], gapaR[i-1][j-1], gapbR[i-1][j-1]) + MISMATCH_SCORE
			}
			// clip of contig until j best mapped to i loc in ref
			if scoreR[i][j] > scoreR[bestR[j]][j] {
				bestR[j] = i
			}
		}
	}

	// backtrace
	max := scoreL[bestL[len(read)]][len(read)]
	split := len(read)
	// find best split
	for i := len(read) - 1; i >= 0; i-- {
		if scoreL[bestL[i]][i]+scoreR[bestR[len(read)-i]][len(read)-i] > max {
			max = scoreL[bestL[i]][i] + scoreR[bestR[len(read)-i]][len(read)-i]
			split = i
		}
	}
	endL := bestL[split] - 1
	result.lbp = endL
	var startR int
	if svtype == "del" {
		startR = len(refR) - bestR[len(read)-split]
		result.rbp = len(refR) - bestR[len(read)-split] - 1
	} else if svtype == "INV" {
		startR = bestR[len(read)-split]
		result.rbp = bestR[len(read)-split] - 1
	}

	aL.Reset()
	bL.Reset()
	cL.Reset()
	cur, match, mismatch, indel := 0, 0, 0, 0

	pj := split        // read iterator
	pi := bestL[split] // ref iterator

	// left part alignment result
	for pi > 0 && pj > 0 {
		if cur == 0 {
			var tmp int
			if refL[pi-1] == read[pj-1] {
				tmp = MATCH_SCORE
				cL.WriteString("|")
			} else {
				tmp = MISMATCH_SCORE
				cL.WriteString(" ")
			}
			aL.WriteString(refL[pi-1 : pi])
			bL.WriteString(read[pj-1 : pj])

			if gapaL[pi-1][pj-1]+tmp == scoreL[pi][pj] {
				indel++
				cur = 1
			} else if gapbL[pi-1][pj-1]+tmp == scoreL[pi][pj] {
				indel++
				cur = 2
			}
			pi--
			pj--

			if tmp > 0 {
				match++
			} else {
				mismatch++
			}
		} else if cur == 1 {
			aL.WriteString(refL[pi-1 : pi])
			bL.WriteString("-")
			cL.WriteString(" ")
			if scoreL[pi-1][pj]+GAP_EXTENSION_SCORE+GAP_OPENING_SCORE == gapaL[pi][pj] {
				cur = 0
			}
			pi--
		} else {
			aL.WriteString("-")
			bL.WriteString(read[pj-1 : pj])
			cL.WriteString(" ")
			if scoreL[pi][pj-1]+GAP_EXTENSION_SCORE+GAP_OPENING_SCORE == gapbL[pi][pj] {
				cur = 0
			}
			pj--
		}
	}
	//for soft clips
	for pj > 0 {
		aL.WriteString("-")
		cL.WriteString("S")
		bL.WriteString(read[pj-1 : pj])
		pj--
	}

	result.aL = Reverse(aL.String())
	result.bL = Reverse(bL.String())
	result.cL = Reverse(cL.String())
	// start loc of mapping in ref
	startL := pi
	result.pos = startL

	// find identity of left part
	l := 1 + math.Abs(float64(endL-startL+1-split))
	result.identityL = -1.0

	if split != 0 {
		logl := math.Log(l)
		result.identityL = 1 - (float64(mismatch+indel)+logl)/float64(split)
	}

	pi = bestR[len(read)-split] // ref iterator
	pj = len(read) - split      // read iterator
	match, mismatch, indel = 0, 0, 0

	aR.Reset()
	bR.Reset()
	cR.Reset()

	// right part alignment result
	for pi > 0 && pj > 0 {
		if cur == 0 {
			var tmp int
			if refR[pi-1] == readReversed[pj-1] {
				tmp = MATCH_SCORE
				cR.WriteString("|")
			} else {
				tmp = MISMATCH_SCORE
				cR.WriteString(" ")
			}
			aR.WriteString(refR[pi-1 : pi])
			bR.WriteString(readReversed[pj-1 : pj])

			if gapaR[pi-1][pj-1]+tmp == scoreR[pi][pj] {
				indel++
				cur = 1
			} else if gapbR[pi-1][pj-1]+tmp == scoreR[pi][pj] {
				indel++
				cur = 2
			}
			pi--
			pj--

			if tmp > 0 {
				match++
			} else {
				mismatch++
			}
		} else if cur == 1 {
			aR.WriteString(refR[pi-1 : pi])
			bR.WriteString("-")
			cR.WriteString(" ")
			if scoreR[pi-1][pj]+GAP_EXTENSION_SCORE+GAP_OPENING_SCORE == gapaR[pi][pj] {
				cur = 0
			}
			pi--
		} else {
			aR.WriteString("-")
			bR.WriteString(readReversed[pj-1 : pj])
			cR.WriteString(" ")
			if scoreR[pi][pj-1]+GAP_EXTENSION_SCORE+GAP_OPENING_SCORE == gapbR[pi][pj] {
				cur = 0
			}
			pj--
		}
	}
	//for soft clips
	for pj > 0 {
		aR.WriteString("-")
		cR.WriteString("S")
		bR.WriteString(readReversed[pj-1 : pj])
		pj--
	}
	endR := len(refR) - pi - 1

	result.aR = aR.String()
	result.bR = bR.String()
	result.cR = cR.String()

	// find identity of right part
	len := len(read) - split
	l = 1 + math.Abs(float64(endR-startR+1-len))
	result.identityR = -1.0
	if len != 0 {
		logl := math.Log(l)
		result.identityR = 1 - (float64(mismatch+indel)+logl)/float64(len)
	}

	return result
}

func max2(a int, b int) int {
	if a > b {
		return a
	}
	return b
}
func max3(a int, b int, c int) int {
	if a > b && a > c {
		return a
	}
	if b > a && b > c {
		return b
	}
	return c
}
