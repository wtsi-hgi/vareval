/* Copyright (c) 2018 Genome Research Ltd.

Authors:
* Sarah Chacko <sc35@sanger.ac.uk>

This file is part of the vareval project.

vareval is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------------------
*/

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
	fmt.Println("Data Checks")
	v := "/home/sjc/Downloads/h38_megaELGH.wg.max_missingness_2pc.with_AF.sorted.vcf.gz"
	r := "/home/sjc/testdata/ROH_regions/allROH.txt"
	m := "/home/sjc/Downloads/sample_id_mappings_egan_to_elgh"
	a, b, c := dataChecks(v, r, m)
	fmt.Println(a, b, c)
}

func testSampleBEDsFromROH(t *testing.T) {
	fmt.Println("SampleBEDsFromROH")
	file1 := "/home/sjc/Downloads/sample_id_mappings_egan_to_elgh"
	rohfile := "/home/sjc/testdata/ROH_regions/allROH.txt"

	m, _, err := mapSampleNames(file1)
	fmt.Println(m)
	if err != nil {
		t.Errorf("%s", err.Error())
	}

	err = sampleBEDsFromROH(rohfile, true, true, m, m)
	if err != nil {
		t.Errorf("%s", err.Error())
	}

}

func TestOverlap(t *testing.T) {
	fmt.Println("Overlap")
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

func TestIndex(t *testing.T) {
	vcf := "use_h38_megaELGH.wg.max_missingness_2pc.with_AF.sorted.vcf.gz"
	err := indexVCF(vcf)
	if err != nil {
		t.Errorf("%s", err.Error())
	}
}
