package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

func main() {
	log.SetFlags(log.LstdFlags | log.Lshortfile)
	// set up parameters
	mapFile := flag.String("map", "", "mapping between sample names")
	//a mapping file of sample names (different in the vcf and the roh files)
	//a list of ROH regions by sample and chromosome
	rohFile := flag.String("roh", "", "file of regions of homozygosity per sample")
	//a multisample vcf
	vcfFile := flag.String("vcf", "", "a multisample vcf to analyse")
	flag.Parse()

	// set up the mapping between sample names (ROH to vcf)
	m, err := mapSampleNames(*mapFile)
	if err != nil {
		log.Println(err)
		return
	}
	// Find the samples from the multisample vcf

	s, err := SamplesFromVCF(*vcfFile)
	if err != nil {
		log.Println(err)
		return
	}
	// make the separate sample bed files from the ROH input file
	err = sampleBEDsFromROH(*rohFile, m)
	if err != nil {
		log.Println(err)
		return
	}

	err = handleSamples(*vcfFile, s)
	if err != nil {
		log.Println(err)
		return
	}
}

// handleSamples works out eight counts for each sample and saves them to a stats file
// format sampleid, allcalls, all het calls,all filtered calls, all filtered het calls,
//  all calls in ROH, all het calls in ROH, all filtered calls in ROH and all filtered het calls in ROH
// if no filter was applied the filtered set will match the unfiltered set.
// the final vcfs of het calls in ROH regions are also saved
func handleSamples(vcf string, s []string) (err error) {
	// base filename
	vcfBase := strings.TrimSuffix(filepath.Base(vcf), ".vcf.gz")
	// open the stats file
	f, err := os.Create("stats_" + vcfBase + ".csv")
	if err != nil {
		return
	}
	defer f.Close()
	//	f.WriteString("sample,allVariants,hetVariants,allFilteredVariants, hetFilteredVariants, allROHVariants,hetROHVariants, allFilteredROHVariants, hetFilteredROHVariants\n")
	f.WriteString("sample,allROHVariants,hetROHVariants,allFilteredROHVariants,hetFilteredROHVariants\n")

	// deal with each sample
	for i := range s {

		// get the single sample vcf
		ssVCF := vcfBase + "_" + s[i] + ".vcf"
		err = singleSampleVCF(s[i], vcf, ssVCF)
		if err != nil {

			return
		}
		// count all calls and het calls
		/*
			all, het, allF, hetF, err := countVariants(ssVCF)
			if err != nil {
				return err
			}
			f.WriteString(fmt.Sprintf("%s, %s,%s,%s,%s", s[i], all, het, allF, hetF))
		*/
		// compress and index
		comSSVCF, err := compressVCF(ssVCF)
		if err != nil {

			return err
		}
		os.Remove(ssVCF)

		_, err = indexVCF(comSSVCF)
		if err != nil {
			return err
		}

		// intersect with roh regions
		intersect, err := intersectVCF(comSSVCF, s[i]+".bed")
		if err != nil {

			return err
		}
		os.Remove(comSSVCF)
		allR, hetR, allRF, hetRF, err := countVariants(intersect)
		if err != nil {
			return err
		}
		f.WriteString(fmt.Sprintf("%s,%s,%s,%s,%s\n", s[i], allR, hetR, allRF, hetRF))

		// output actual het calls while testing
		err = hetCalls(intersect, "het_"+intersect)

		if err != nil {
			return err
		}
		os.Remove(intersect)

	}
	return
}

// SamplesFromVCF will error if can't find bcftools .. will be run in docker containing it
func SamplesFromVCF(vcf string) (samples []string, err error) {
	//vcf = "/home/sjc/Downloads/ELGH-V2plus.gvcf_to_vcf.20180516.vcf.gz"
	s, err := exec.Command("bcftools", "query", "-l", vcf).Output()
	//err := cmd1.Run()
	//log.Printf("Command finished with output %v and error: %v", string(s), err)
	if err != nil {
		err = fmt.Errorf("Couldn't get samples from vcf, is bcftools installed? %s", err.Error())
		return
	}
	samples = strings.Fields(string(s))
	return
}

func singleSampleVCF(sample string, vcf string, outputvcf string) (err error) {
	//vcfBase := strings.TrimSuffix(filepath.Base(vcf), ".vcf.gz")
	//file = vcfBase + "_" + sample + ".vcf"

	_, err = exec.Command("bcftools", "view", "-s", sample, "-c1", "-o", outputvcf, vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't get calls for sample %s from vcf %s, is bcftools installed? %s", sample, vcf, err.Error())
		return
	}
	return
}

func compressVCF(vcf string) (file string, err error) {
	file = vcf + ".gz"

	_, err = exec.Command("bcftools", "view", "-Oz", "-o", file, vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't compress  vcf %s, is bcftools installed? %s", vcf, err.Error())
		return
	}
	return
}

func indexVCF(vcf string) (file string, err error) {
	file = vcf + ".gz"

	_, err = exec.Command("bcftools", "index", vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't index  vcf %s, is bcftools installed? %s", vcf, err.Error())
		return
	}
	return
}

func intersectVCF(vcf string, bed string) (file string, err error) {
	file = "in_roh_" + strings.TrimSuffix(vcf, ".gz")

	_, err = exec.Command("bcftools", "view", "-R", bed, "-o", file, vcf).Output()

	if err != nil {
		err = fmt.Errorf("Couldn't intersect  vcf %s with bed %s is bcftools installed? %s", vcf, bed, err.Error())
		return
	}

	return
}

func hetCalls(vcf string, outvcf string) (err error) {
	//get het calls with filters applied for a single sample

	_, err = exec.Command("bcftools", "view", "-g", "het", "-o", outvcf, vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't get het calls for file %s, is bcftools installed? %s", vcf, err.Error())
		return
	}
	return
}

func mapSampleNames(mapfile string) (m map[string]string, err error) {
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
		if len(next2) > 1 {
			m[next2[1]] = next2[0]
		}
	}
	return
}

// count all variants, het variants, all filtered variants and het filtered variants from
// a vcf (there may be no filter so use  -f .,PASS)
// apply the filters or intersect, take off the header, count the lines
func countVariants(vcf string) (all string, het string, allFiltered string, hetsFiltered string, err error) {
	headerlessVCF := "all_headless_" + vcf
	// all, unfiltered
	_, err = exec.Command("bcftools", "view", "-H", "-o", headerlessVCF, vcf).Output()
	if err != nil {
		return
	}
	out, err := exec.Command("wc", "-l", headerlessVCF).Output()
	if err != nil {
		return
	}
	all = strings.Split(string(out), " ")[0]

	// hets, unfiltered
	allHets := "hets" + vcf
	headerlessAllHets := "hets_headless" + allHets
	err = hetCalls(vcf, allHets)
	if err != nil {

		return
	}
	_, err = exec.Command("bcftools", "view", "-H", "-o", headerlessAllHets, allHets).Output()
	if err != nil {
		return
	}

	out, err = exec.Command("wc", "-l", headerlessAllHets).Output()
	if err != nil {
		return
	}

	het = strings.Split(string(out), " ")[0]

	// all, filtered
	_, err = exec.Command("bcftools", "view", "-H", "-f", ".,PASS", "-o", headerlessVCF, vcf).Output()
	if err != nil {
		return
	}
	out, err = exec.Command("wc", "-l", headerlessVCF).Output()
	if err != nil {
		return
	}

	allFiltered = strings.Split(string(out), " ")[0]

	// hets filtered
	_, err = exec.Command("bcftools", "view", "-H", "-f", ".,PASS", "-o", "filtered"+allHets, allHets).Output()
	if err != nil {
		return
	}

	out, err = exec.Command("wc", "-l", "filtered"+allHets).Output()
	if err != nil {
		return
	}

	hetsFiltered = strings.Split(string(out), " ")[0]

	os.Remove(headerlessVCF)
	os.Remove(allHets)
	os.Remove(headerlessAllHets)
	os.Remove("filtered" + allHets)

	if err != nil {
		fmt.Println(err)
		return
	}

	return

}

func sampleBEDsFromROH(roh string, m map[string]string) (err error) {
	// the ROH has format RG, sample, chromosome, start, end, length, numberofMarkers, Quality
	// the sample does not match the sample in the vcf, the mapping file is needed.
	// the chromosome may not be the same as in the vcf, can be chr1 or 1 etc

	// we need, for a sample, a bed file of chromosome, start, end
	bed := make(map[string][]string)
	f, err := os.Open(roh)
	if err != nil {
		err = fmt.Errorf("Couldn't open the roh file: %s", err.Error())
		return
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		l := strings.TrimSpace(scanner.Text())
		if strings.HasPrefix(l, "#") {
			continue // header
		}
		line := strings.Fields(l)
		if len(line) < 5 {
			continue // not a real line
		}
		chromosome := line[2]
		if !strings.HasPrefix(chromosome, "chr") {
			//fmt.Println("Converted chromosome names in ROH file")
			chromosome = "chr" + chromosome
		}
		sample1 := line[1]
		sample := ""
		if val, ok := m[sample1]; ok {
			sample = val
		} else {
			sample = "unmapped_" + sample1
			//fmt.Println(fmt.Sprintf("sample in regions not in vcf %s (OK)", sample1))
		}
		if strings.TrimSpace(sample) == "" {
			sample = "blank_" + sample1
		}
		start := line[3]
		end := line[4]
		if (!strings.HasPrefix(sample, "unmapped")) && (!strings.HasPrefix(sample, "blank")) {
			// only add samples from the vcf to the bed files
			bedfields := []string{chromosome, start, end, sample, sample1} /// samples ignored in processing
			bedline := strings.Join(bedfields, "\t")
			//fmt.Println(bedline)

			if val, ok := bed[sample]; ok {
				newval := append(val, bedline)
				bed[sample] = newval
			} else {
				bed[sample] = []string{bedline}
			}
		}

	}

	for k, v := range bed {
		f, err := os.Create(k + ".bed")
		if err != nil {
			return err
		}
		defer f.Close()
		sort.Sort(bedLine(v))
		for i := range v {
			f.WriteString(v[i] + "\n")
		}
	}

	return
}

// sorting util, sort bed file by chromosome (field 1) and by start position
type bedLine []string

func (s bedLine) Len() int {
	return len(s)
}
func (s bedLine) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}
func (s bedLine) Less(i, j int) bool {
	fi := strings.Fields(strings.TrimSpace(string(s[i])))
	fj := strings.Fields(strings.TrimSpace(string(s[j])))
	chromosomei := fi[0]
	starti := fi[1]

	if strings.HasPrefix(chromosomei, "chr") {
		chromosomei = chromosomei[3:]
	}
	chromosomej := fj[0]
	startj := fj[1]
	if strings.HasPrefix(chromosomej, "chr") {
		chromosomej = chromosomej[3:]
	}

	var firstIsLess bool
	ci, _ := strconv.Atoi(chromosomei)
	cj, _ := strconv.Atoi(chromosomej)

	switch {
	case ci < cj:
		firstIsLess = true
	case ci == cj:
		si, _ := strconv.Atoi(starti)
		sj, _ := strconv.Atoi(startj)
		if si < sj {
			firstIsLess = true
		} else {
			firstIsLess = false
		}
	case i > cj:
		firstIsLess = false
	}

	return firstIsLess
}
