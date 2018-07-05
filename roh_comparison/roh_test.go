package main

import (
	"fmt"
	"testing"
)

func testSamplesAndChrFromROH(t *testing.T) {
	v := "/home/sjc/testdata/ROH_regions/allROH.txt"
	a, b, c := SamplesAndChrFromROH(v)
	fmt.Println(a, b, c)

}

func testSamplesAndChrFromVCF(t *testing.T) {
	v := "/home/sjc/Downloads/h38_megaELGH.wg.max_missingness_2pc.with_AF.sorted.vcf.gz"
	a, b, c := SamplesAndChrFromVCF(v)
	fmt.Println(a, b, c)
}

func TestDataChecks(t *testing.T) {
	v := "/home/sjc/Downloads/h38_megaELGH.wg.max_missingness_2pc.with_AF.sorted.vcf.gz"
	r := "/home/sjc/testdata/ROH_regions/allROH.txt"
	m := "/home/sjc/Downloads/sample_id_mappings_egan_to_elgh"
	a, b, c := dataChecks(v, r, m)
	fmt.Println(a, b, c)
}

func testSampleBEDsFromROH(t *testing.T) {
	file1 := "/home/sjc/Downloads/sample_id_mappings_egan_to_elgh"
	rohfile := "/home/sjc/testdata/ROH_regions/allROH.txt"

	m, _, err := mapSampleNames(file1)
	fmt.Println(m)
	if err != nil {
		t.Errorf("%s", err.Error())
	}

	err = sampleBEDsFromROH(rohfile, m)
	if err != nil {
		t.Errorf("%s", err.Error())
	}

}

func TestOverlap(t *testing.T) {
	m := make(map[string]bool)
	m["a"] = true
	m["b"] = true
	m["c"] = true
	m["d"] = true
	m["e"] = true
	n := make(map[string]bool)
	/*n["a"] = true
	n["b"] = true
	n["d"] = true
	n["x"] = true*/
	fmt.Println(overlap(m, n))
}
