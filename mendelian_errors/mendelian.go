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

// This package uses known 'truth' data about relatedness (mother/father/child trios) in samples
// to find variants that could not be inherited (using bcftools mendelian plugin)
// so are either errors or de novo mutations,
// as a check on the accuracy of different variant callers
package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"
)

// ./mendelian_errors -map= -trios= -vcf=
func main() {
	log.SetFlags(log.LstdFlags | log.Lshortfile)
	// set up parameters
	// the defaults are for ease of use when run in docker, you map the real files
	// to the default file names.
	mapFile := flag.String("map", "/tmp/map.txt", "mapping between stable ids and sample names")
	//a mapping file of sample names (different in the vcf and the trios files)
	//a list of ROH regions by sample and chromosome
	triosFile := flag.String("trios", "/tmp/alltrios.txt", "file of trio relationships (using stable ids)")
	//a multisample vcf
	vcfFile := flag.String("vcf", "/tmp/vcf.vcf.gz", "a multisample vcf to analyse")
	flag.Parse()

	triosTempFile := "vcf_trios.txt"

	// Find the samples from the multisample vcf
	s, err := SamplesFromVCF(*vcfFile)
	if err != nil {
		log.Println(err)
		os.Exit(2)
	}

	// set up the mapping between sample names (trios to vcf)
	m, err := mapSampleNames(*mapFile, s)
	if err != nil {
		log.Println(err)
		os.Exit(1)
	}

	// make the corrected trios file for these samples
	err = makeTriosFile(*triosFile, triosTempFile, m, s)
	if err != nil {
		log.Println(err)
		os.Exit(3)
	}

	err = mendelianPlugin(*vcfFile, triosTempFile)
	if err != nil {
		log.Println(err)
		os.Exit(4)
	}

}

func mendelianPlugin(vcf string, trios string) (err error) {
	// call the bcf plugin
	// bcftools +mendelian <vcf>  -T <trios> -c > <output>
	// the -- changed between 1.3.1 and 1.7 not sure where
	cmd := exec.Command("bcftools", "+mendelian", vcf /*"--",*/, "-T", trios, "-c")
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	err = cmd.Run()
	if err != nil {
		log.Fatalf("Failed to call bcftools mendelian plugin: cmd.Run() failed with %s\n", err)
	}
	outStr, errStr := string(stdout.Bytes()), string(stderr.Bytes())
	fmt.Printf("out:\n%s\nerr:\n%s\n", outStr, errStr)
	/* changed as the output is in stderr and this is not a good way to capture it
	output, err := exec.Command("bcftools", "+mendelian", vcf, "-T", trios, "-c", "2>xx").Output()
	if err != nil {
		err = fmt.Errorf("Couldn't run mendelian plugin for vcf %s with trios %s : %s", vcf, trios, err.Error())
		return
	}
	fmt.Println(string(output))*/
	f, err := os.Create("mendelian_output")
	if err != nil {
		err = fmt.Errorf("Couldn't open the output file : %s", err.Error())
		return
	}
	defer f.Close()

	f.WriteString(string(errStr))

	return
}

// SamplesFromVCF will error if can't find bcftools .. will be run in docker containing it
// using a map as we need to search it, easier
func SamplesFromVCF(vcf string) (samples map[string]bool, err error) {

	samples = make(map[string]bool) //vcf = "/home/sjc/Downloads/ELGH-V2plus.gvcf_to_vcf.20180516.vcf.gz"
	s, err := exec.Command("bcftools", "query", "-l", vcf).Output()
	//err := cmd1.Run()
	//log.Printf("Command finished with output %v and error: %v", string(s), err)
	if err != nil {
		err = fmt.Errorf("Couldn't get samples from vcf, is bcftools installed? %s", err.Error())
		return
	}
	temp := strings.Fields(string(s))
	for i := range temp {
		samples[temp[i]] = true
	}
	return
}

// mapSampleNames sets up a mapping from stable ids, used in original trios file
// to ega ids used in the vcf. Each line of the mapping file has
// decipher_id	person_stable_id	sanger_id	ega_id	is_proband	gender
// and we want to map person_stable_id to ega_id
// ??? could use gender as a check?
func mapSampleNames(mapfile string, samples map[string]bool) (m map[string]string, err error) {
	m = make(map[string]string)
	f, err := os.Open(mapfile)
	if err != nil {
		err = fmt.Errorf("Couldn't open the mapping file for sample names: %s", err.Error())
		return
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		next := strings.TrimSpace(scanner.Text())
		next2 := strings.Fields(next)
		if len(next2) >= 6 {
			if _, ok := samples[next2[3]]; ok { // only map samples we need
				m[next2[1]] = next2[3]
			}
		} else {
			err = fmt.Errorf("Unexpected format for mapping file %s, should have 6 fields", mapfile)
			return
		}
	}
	return
}

// make a trios file (father mother child) for the sample ids used in the vcf
// with the trios data from the input trios file (another id then the three required )
func makeTriosFile(inFile string, outfile string, m map[string]string, samples map[string]bool) (err error) {
	fOut, err := os.Create(outfile)
	if err != nil {
		err = fmt.Errorf("Couldn't create the new trios file: %s", outfile)
		return
	}
	defer fOut.Close()
	fIn, err := os.Open(inFile)
	if err != nil {
		err = fmt.Errorf("Couldn't open the existing trios file: %s", inFile)
		return
	}
	defer fIn.Close()

	scanner := bufio.NewScanner(fIn)

	for scanner.Scan() {
		next := strings.TrimSpace(scanner.Text())
		next2 := strings.Fields(next)
		if len(next2) != 4 {
			err = fmt.Errorf("Unexpected format for trios file %s, should have 4 fields", inFile)
			return
		}

		nextLine := fmt.Sprintf("%s,%s,%s\n", m[next2[1]], m[next2[2]], m[next2[3]])
		_, ok1 := samples[m[next2[1]]]
		_, ok2 := samples[m[next2[2]]]
		_, ok3 := samples[m[next2[3]]]
		if ok1 && ok2 && ok3 { // all of trio exists in our dataset
			fOut.WriteString(nextLine)
		}
	}
	return
}
