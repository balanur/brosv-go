package main

type SV struct {
	id         string
	Chromosome string
	Start      int
	End        int
	Type       string
}

type SVStore struct {
	svMap map[string]SV
}

func NewSVStore() SVStore {
	var result SVStore
	result.svMap = make(map[string]SV)
	return result
}

func (svStore *SVStore) add(sv SV) {
	svStore.svMap[sv.id] = sv
}

func (svStore *SVStore) get(id string) SV {
	return svStore.svMap[id]
}

type Interval struct {
	head int
	tail int
	svId string
	side Side
}

type CIStore struct {
	ciList []Interval
	ciMap  map[string][]int
}

func NewCIStore() CIStore {
	var result CIStore
	result.ciMap = make(map[string][]int)
	return result
}

func (ciStore *CIStore) add(svStore SVStore, interval Interval) {
	ciStore.ciList = append(ciStore.ciList, interval)
	ciIndex := len(ciStore.ciList) - 1
	chrName := svStore.get(interval.svId).Chromosome
	ciStore.ciMap[chrName] = append(ciStore.ciMap[chrName], ciIndex)
}

func (ciStore *CIStore) get(index int) Interval {
	return ciStore.ciList[index]
}

func (ciStore *CIStore) findIntersectingIntervals(chrName string, start int, end int) []int {
	intervalsInChromosome := ciStore.ciMap[chrName]
	var result []int
	for _, i := range intervalsInChromosome {
		interval := ciStore.get(i)
		if (interval.tail >= start && start >= interval.head) || (interval.tail >= end && end >= interval.head) {
			result = append(result, i)
		}
	}
	return result
}

type IndexPair struct {
	mappingOri int
	pairNumber int
	ciIndex    int
}

const (
	RP_CONC = iota + 1
	RP_UNMAPPED
	RP_MAPPED
	RP_BOTH
)

// Side of Confidence interval
type Side int

const (
	leftCI Side = iota + 1
	rightCI
)

// FaiEntry is a line from fasta index
type FaiEntry struct {
	title     string
	length    int64
	offset    int64
	linebases int
	linewidth int
}

// Chromosome title and content
type Chromosome struct {
	title   string
	content string
}

// Genome : Whole genome consisting of fasta index and chromosome contents
type Genome struct {
	faiEntries []FaiEntry
	chms       []Chromosome
	chmMap     map[string]int
}

func (genome *Genome) getChr(title string) Chromosome {
	return genome.chms[genome.chmMap[title]]
}
