package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"
)

func main() {
	log.SetFlags(log.LstdFlags | log.Lshortfile)
	// set up parameters
	mapFile := flag.String("map", "", "mapping between stable ids and sample names")
	//a mapping file of sample names (different in the vcf and the trios files)
	//a list of ROH regions by sample and chromosome
	triosFile := flag.String("trios", "", "file of trio relationships (using stable ids)")
	//a multisample vcf
	vcfFile := flag.String("vcf", "", "a multisample vcf to analyse")
	flag.Parse()

	// Find the samples from the multisample vcf
	s, err := SamplesFromVCF(*vcfFile)
	if err != nil {
		log.Println(err)
		os.Exit(2)
	}

	// set up the mapping between sample names (ROH to vcf)
	m, err := mapSampleNames(*mapFile, s)
	if err != nil {
		log.Println(err)
		os.Exit(1)
	}

	// make the corrected trios file for these samples
	err = makeTriosFile(*triosFile, "vcf_trios.txt", m, s)
	if err != nil {
		log.Println(err)
		os.Exit(3)
	}

	err = mendelianPlugin(*vcfFile, "vcf_trios.txt")
	if err != nil {
		log.Println(err)
		os.Exit(4)
	}

}

func mendelianPlugin(vcf string, trios string) (err error) {
	// call the bcf plugin
	// bcftools +mendelian <vcf>  -T <trios> -c > <output>
	output, err := exec.Command("bcftools", "+mendelian", vcf, "-T", trios, "-c").Output()
	if err != nil {
		log.Println(err)
		return
	}
	fmt.Println(string(output))
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

		nextLine := fmt.Sprintf("%s\t%s\t%s\n", m[next2[1]], m[next2[2]], m[next2[3]])
		_, ok1 := samples[m[next2[1]]]
		_, ok2 := samples[m[next2[2]]]
		_, ok3 := samples[m[next2[3]]]
		if ok1 && ok2 && ok3 { // all of trio exists in our dataset
			fOut.WriteString(nextLine)
		}
	}
	return
}
