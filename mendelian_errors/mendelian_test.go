/* Copyright (c) 2018 Genome Research Ltd.

Authors:
* Sarah Chacko sc35@sanger.ac.uk

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
	"os"
	"testing"
)

func TestMapSampleNames(t *testing.T) {
	// decipher_id	person_stable_id	sanger_id	ega_id	is_proband	gender
	temp, err := os.Create("temp.txt")
	if err != nil {
		t.Errorf("%s", err.Error())
	}

	temp.WriteString(fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", "a", "b", "c", "d", "e", "f"))
	temp.WriteString(fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", "aa", "ba", "ca", "da", "ea", "fa"))
	temp.Close()
	samples := make(map[string]bool)
	samples["d"] = true
	samples["da"] = true
	samples["ccc"] = true
	samples["dcc"] = true

	m, err := mapSampleNames("temp.txt", samples)
	if m["b"] != "d" {
		t.Errorf("Expected %s got %s", "d", m["b"])
	}
	if len(m) != 2 {
		t.Errorf("Expected %d got %d", 2, len(m))
	}

	os.Remove("temp.txt")
}

func TestMakeTriosFile(t *testing.T) {
	// decipher_id	person_stable_id	sanger_id	ega_id	is_proband	gender
	temp, err := os.Create("temp.txt")
	if err != nil {
		t.Errorf("%s", err.Error())
	}

	temp.WriteString(fmt.Sprintf("%s\t%s\t%s\n", "a", "b", "c"))
	temp.WriteString(fmt.Sprintf("%s\t%s\t%s\n", "da", "ea", "fa"))
	temp.Close()

	m := make(map[string]string)
	m["a"] = "z"
	m["b"] = "y"
	m["c"] = "x"
	m["da"] = "xz"
	m["ea"] = "zz"
	m["fa"] = "fz"

	s := make(map[string]bool)
	s["x"] = true
	s["y"] = true
	s["z"] = true
	s["xz"] = true
	s["zz"] = true
	s["fz"] = true

	makeTriosFile("temp.txt", "tr.txt", m, s)

}

func TestMendelianPlugin(t *testing.T) {
	// takes a multisample vcf and a trios file and outputs counts of mendelian errors
}
