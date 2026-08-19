// Native baseline for the CJR zero-vertex classification hot loop.
package main

import (
	"flag"
	"fmt"
	"time"
)

type caseData struct {
	genera       []int
	degrees      []int
	edgeZeros    []int
	markingZeros []int
}

func cases(count int) []caseData {
	state := uint32(0xC0FFEE)
	draw := func() uint32 {
		state = 1664525*state + 1013904223
		return state
	}
	answer := make([]caseData, 0, count)
	for index := 0; index < count; index++ {
		zeroCount := 1 + int(draw()%8)
		item := caseData{
			genera:  make([]int, zeroCount),
			degrees: make([]int, zeroCount),
		}
		for zero := 0; zero < zeroCount; zero++ {
			if draw()%17 == 0 {
				item.genera[zero] = 1
			}
		}
		for zero := 0; zero < zeroCount; zero++ {
			if draw()%19 == 0 {
				item.degrees[zero] = 1
			}
		}
		edgeCount := 1 + int(draw()%12)
		item.edgeZeros = make([]int, edgeCount)
		for edge := range item.edgeZeros {
			item.edgeZeros[edge] = int(draw() % uint32(zeroCount))
		}
		markingCount := int(draw() % 5)
		item.markingZeros = make([]int, markingCount)
		for marking := range item.markingZeros {
			item.markingZeros[marking] = int(draw() % uint32(zeroCount))
		}
		answer = append(answer, item)
	}
	return answer
}

func classify(item caseData) int {
	zeroCount := len(item.genera)
	edgeValences := make([]int, zeroCount)
	markingValences := make([]int, zeroCount)
	for _, zero := range item.edgeZeros {
		edgeValences[zero]++
	}
	for _, zero := range item.markingZeros {
		markingValences[zero]++
	}
	checksum := 0
	for zero := 0; zero < zeroCount; zero++ {
		edgeValence := edgeValences[zero]
		markingValence := markingValences[zero]
		valence := edgeValence + markingValence
		kind := 5
		if item.genera[zero] != 0 || item.degrees[zero] != 0 || valence > 2 {
			kind = 1
		} else if edgeValence == 1 && markingValence == 0 {
			kind = 2
		} else if edgeValence == 1 && markingValence == 1 {
			kind = 3
		} else if edgeValence == 2 && markingValence == 0 {
			kind = 4
		}
		checksum += kind * (zero + 1)
	}
	return checksum
}

func main() {
	caseCount := flag.Int("cases", 4096, "number of deterministic cases")
	rounds := flag.Int("rounds", 500, "number of passes over the cases")
	flag.Parse()
	workload := cases(*caseCount)
	started := time.Now()
	checksum := 0
	for round := 0; round < *rounds; round++ {
		for _, item := range workload {
			checksum += classify(item)
		}
	}
	elapsed := time.Since(started)
	calls := *caseCount * *rounds
	fmt.Printf("language=go calls=%d seconds=%.6f ns_per_call=%.1f checksum=%d\n",
		calls, elapsed.Seconds(), float64(elapsed.Nanoseconds())/float64(calls), checksum)
}
