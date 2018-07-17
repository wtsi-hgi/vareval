package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"time"
)

func main() {

	// set up parameters
	mapFile := flag.String("map", "", "mapping between sample names")
	//a mapping file of sample names (different in the vcf and the roh files)
	//a list of ROH regions by sample and chromosome
	rohFile := flag.String("roh", "", "file of regions of homozygosity per sample")
	//a multisample vcf
	vcfFile := flag.String("vcf", "", "a multisample vcf to analyse")
	keepHetVCFs := flag.Bool("keep", true, "keep the het call vcfs (as well as stats)")
	flag.Parse()

	// set up logging
	f, err := os.Create("log_" + strings.TrimSuffix(filepath.Base(*vcfFile), ".vcf.gz") + ".txt")
	if err != nil {
		log.Fatalf("error opening file: %v", err)
	}
	defer f.Close()
	log.SetOutput(f)
	log.SetFlags(log.LstdFlags | log.Lshortfile)
	log.Println(fmt.Sprintf("Starting with file %s, mapping file %s and ROH file %s", *vcfFile, *mapFile, *rohFile))

	// check data ... do sample names and chromosome names match in vcf and roh files
	mapSamples, mapChr, err := dataChecks(*vcfFile, *rohFile, *mapFile)
	if err != nil {
		log.Println(err)
		os.Exit(1)
	}
	if mapSamples {
		log.Println("samples in ROH file and vcf do not match, mapping file will be used")
	} else {
		log.Println("sample names match in ROH file and vcf")
	}
	if mapChr {
		log.Println("Chromosome names in ROH file and vcf do not match, chr version will be used")
	} else {
		log.Println("Chromosome names match in ROH file and vcf")
	}

	// set up the mapping between sample names (ROH to vcf)
	var m1, m2 map[string]string
	//if mapSamples { needed anyway for filtering
	m1, m2, err = mapSampleNames(*mapFile)
	if err != nil {
		log.Println(err)
		os.Exit(2)
	}
	log.Println("Set up sample maps")

	//}
	// Find the samples from the multisample vcf which
	// are also in the map file
	s, err := SamplesFromVCF(*vcfFile, m1, m2)
	if err != nil {
		log.Println(err)
		os.Exit(3)
	}
	log.Println(fmt.Sprintf("Got %d samples from vcf", len(s)))
	// make the separate sample bed files from the ROH input file
	err = sampleBEDsFromROH(*rohFile, mapSamples, mapChr, m1, m2)
	if err != nil {
		log.Println(err)
		os.Exit(4)
	}
	log.Println("Sample bed files created")

	//make temp file from wanted samples
	vcf := "use_" + filepath.Base(*vcfFile)
	err = multiSampleVCF(s, *vcfFile, "use_"+filepath.Base(*vcfFile))
	if err != nil {
		log.Println(err)
		//os.Exit(5)
		// if this fails, just use the input vcf
		vcf = *vcfFile
	}
	log.Println("Subset of vcf created with required samples")

	err = indexVCF(vcf)
	if err != nil {
		log.Println(err)
		//os.Exit(6)// usual error is index already there 255
	}
	log.Println("Subset of vcf indexed")

	// handle each sample
	log.Println("Starting sample handling")
	err = handleSamples(vcf, s, m1, m2, *keepHetVCFs)
	if err != nil {
		log.Println(err)
		os.Exit(7)
	}

	log.Println("Finished")
}

// handleSamples works out eight counts for each sample and saves them to a stats file
// format sampleid, allcalls, all het calls,all filtered calls, all filtered het calls,
//  all calls in ROH, all het calls in ROH, all filtered calls in ROH and all filtered het calls in ROH
// if no filter was applied the filtered set will match the unfiltered set.
// the final vcfs of het calls in ROH regions are also saved
func handleSamples(vcf string, s []string, m1, m2 map[string]string, keep bool) (err error) {
	// base filename
	vcfBase := strings.TrimSuffix(vcf, ".vcf.gz")
	// open the stats file
	f, err := os.Create("stats_" + vcfBase + ".csv")
	if err != nil {
		return
	}
	defer f.Close()
	f.WriteString("sample, sample2, hetFilteredVariants, allROHVariants,hetROHVariants, allFilteredROHVariants, hetFilteredROHVariants\n")
	//f.WriteString("sample,allROHVariants,hetROHVariants,allFilteredROHVariants,hetFilteredROHVariants\n")

	// deal with each sample
	for i := range s {
		fmt.Printf("%s: Sample %s: %d of %d\n", time.Now().String(), s[i], i, len(s))
		if _, err := os.Stat(s[i] + ".bed"); os.IsNotExist(err) {
			// does not exist in ROH so no bed file- next sample
			fmt.Println("no bed file for", s[i])
			continue
		}

		// get the single sample vcf
		ssVCF := vcfBase + "_" + s[i] + ".vcf.gz"
		err = singleSampleVCF(s[i], vcf, ssVCF)
		if err != nil {

			return
		}
		// count all calls and het calls slow if used
		/*
			all, het, allF, hetF, err := countVariants(ssVCF)
			if err != nil {
				return err
			}
			f.WriteString(fmt.Sprintf("%s, %s,%s,%s,%s", s[i], all, het, allF, hetF))
		*/
		// compress and index
		/*comSSVCF, err := compressVCF(ssVCF)
		if err != nil {

			return err
		}*/
		comSSVCF := ssVCF // temp. was making then compressing
		//os.Remove(ssVCF)

		err = indexVCF(comSSVCF)
		if err != nil {
			fmt.Println(err)
			// but carry on ... usually this fails with the index existing
		}

		// intersect with roh regions if we have this sample
		if _, err := os.Stat(s[i] + ".bed"); os.IsNotExist(err) {
			// does not exist - not a sample we are using
			continue
		}
		intersect, err := intersectVCF(comSSVCF, s[i]+".bed")
		if err != nil {

			return err
		}
		os.Remove(comSSVCF)
		os.Remove(comSSVCF + ".csi")
		os.Remove(comSSVCF + ".tbi")
		allR, hetR, allRF, hetRF, err := countVariants(intersect)
		if err != nil {
			return err
		}
		// both sample names
		s1 := s[i]
		s2 := s[i]
		if v, ok := m1[s[i]]; ok {
			s1 = v
		}
		if v, ok := m2[s[i]]; ok {
			s2 = v
		}
		f.WriteString(fmt.Sprintf("%s,%s, %s,%s,%s,%s\n", s1, s2, allR, hetR, allRF, hetRF))
		// replace by following line if also finding all calls (ie not in ROH regions)
		//f.WriteString(fmt.Sprintf(",%s,%s,%s,%s\n", allR, hetR, allRF, hetRF))

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
// only use samples we have in the mapping file but report others
func SamplesFromVCF(vcf string, m1 map[string]string, m2 map[string]string) (samples []string, err error) {
	s, err := exec.Command("bcftools", "query", "-l", vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't get samples from vcf, is bcftools installed? %s", err.Error())
		return
	}
	s2 := strings.Fields(string(s))
	// limit to samples we include in the map file
	// but could be either name
	for i := range s2 {
		if _, ok := m1[s2[i]]; ok { /// exists in forward map
			samples = append(samples, s2[i])
		} else if _, ok := m2[s2[i]]; ok { /// exists in freverse map
			samples = append(samples, s2[i])
		}

	}
	s = nil
	s2 = nil
	runtime.GC() // see whether this helps with slowing down

	return
}

// subset vcf by sample list
func multiSampleVCF(samples []string, vcf string, outputvcf string) (err error) {
	//vcfBase := strings.TrimSuffix(filepath.Base(vcf), ".vcf.gz")
	//file = vcfBase + "_" + sample + ".vcf"

	_, err = exec.Command("bcftools", "view", "-s", strings.Join(samples, ","), "-c1", "-o", outputvcf, "-O", "z", vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't get calls for samples %s from vcf %s, is bcftools installed? %s", strings.Join(samples, ","), vcf, err.Error())

	}
	return
}

func singleSampleVCF(sample string, vcf string, outputvcf string) (err error) {
	//vcfBase := strings.TrimSuffix(filepath.Base(vcf), ".vcf.gz")
	//file = vcfBase + "_" + sample + ".vcf"

	_, err = exec.Command("bcftools", "view", "-s", sample, "-c1", "-o", outputvcf, "-Oz", vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't get calls for sample %s from vcf %s, is bcftools installed? %s", sample, vcf, err.Error())

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

func indexVCF(vcf string) (err error) {
	//file = vcf + ".gz"
	// use -f or it fails if index already found
	out, err := exec.Command("bcftools", "index", "-f", vcf).Output()
	if err != nil {
		fmt.Println(out)
		fmt.Println(err)
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

func mapSampleNames(mapfile string) (mForward map[string]string, mReverse map[string]string, err error) {
	mForward = make(map[string]string)
	mReverse = make(map[string]string)
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
			mForward[next2[1]] = next2[0]
			mForward[next2[0]] = next2[1]
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

func sampleBEDsFromROH(roh string, renameSamples, renameChromosomes bool, m1 map[string]string, m2 map[string]string) (err error) {
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
		if renameChromosomes && !strings.HasPrefix(chromosome, "chr") {
			//fmt.Println("Converted chromosome names in ROH file")
			chromosome = "chr" + chromosome
		}
		sample1 := line[1]
		sample := ""
		altSample := ""
		if renameSamples {
			if val, ok := m1[sample1]; ok {

				altSample = sample
				sample = val
			} else {
				sample = "unmapped_" + sample1
				//fmt.Println(fmt.Sprintf("sample in regions not in vcf %s (OK)", sample1))
			}

			if strings.TrimSpace(sample) == "" {
				sample = "blank_" + sample1
			}
		} else { // no mapping but must be in the map file
			if v, ok := m1[sample1]; ok {
				sample = sample1
				altSample = v
			} else if v, ok := m2[sample1]; ok {
				sample = sample1
				altSample = v
				//fmt.Println(fmt.Sprintf("sample in regions not in vcf %s (OK)", sample1))
			} else {
				continue // don't need this one
			}
			//sample = sample1
		}
		start := line[3]
		end := line[4]
		if (!strings.HasPrefix(sample, "unmapped")) && (!strings.HasPrefix(sample, "blank")) {
			// only add real samples from the vcf to the bed files
			bedfields := []string{chromosome, start, end, sample, altSample} /// samples ignored in processing
			bedline := strings.Join(bedfields, "\t")
			//fmt.Println(bedline)
			// for now make both sample sets
			if val, ok := bed[sample]; ok {
				newval := append(val, bedline)
				bed[sample] = newval
				//bed[sample1] = newval
			} else {
				bed[sample] = []string{bedline}
				//bed[sample1] = []string{bedline}
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

func SamplesAndChrFromROH(roh string) (s map[string]bool, chr map[string]bool, err error) {
	s = make(map[string]bool)
	chr = make(map[string]bool)

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
		chr[chromosome] = true
		sample := line[1]
		s[sample] = true

	}
	return
}

func SamplesAndChrFromVCF(vcf string) (s map[string]bool, chr string, err error) {
	s = make(map[string]bool)
	//bcftools query -l x.vcf.gz
	out, err := exec.Command("bcftools", "query", "-l", vcf).Output()
	if err != nil {
		err = fmt.Errorf("Couldn't get samples from vcf, is bcftools installed? %s", err.Error())
		return
	}
	s2 := strings.Fields(string(out))
	for i := range s2 {
		s[s2[i]] = true
	}
	chr = ""
	//bcftools query -f '%CHROM' x.vcf.gz .. too long
	// try to get chromosome names from contig lines in header, else first field of first line
	//bcftools view -H h38_megaELGH.wg.max_missingness_2pc.with_AF.sorted.vcf.gz | head -n 1

	out, err = exec.Command("bcftools", "view", "-H", vcf, "|", "-n", "1").Output()
	if err != nil {
		err = fmt.Errorf("Couldn't get chromosomes from vcf, is bcftools installed? %s", err.Error())
		return
	}
	s2 = strings.Fields(string(out))
	//fmt.Println(s2)
	if strings.HasPrefix("s2[0]", "chr") {
		chr = "chr"
	} else {
		chr = "num"
	}
	//

	return
}

// The sample names should match or need mapping
// the chromosome names should match or need mapping
// all roh samplenames should be in the mapping file if mapping is needed
// report different counts of samples in roh, vcf or mapping file
//
func dataChecks(vcfFile, rohFile, mapFile string) (mapSamples, mapChrs bool, err error) {

	vs, vc, err := SamplesAndChrFromVCF(vcfFile)
	if err != nil {
		err = fmt.Errorf("Datachecks 1: %s", err.Error())
		return
	}
	//fmt.Println(vs, vc)
	rs, rc, err := SamplesAndChrFromROH(rohFile)
	if err != nil {
		err = fmt.Errorf("Datachecks 2: %s", err.Error())
		return
	}
	//	fmt.Println(rs, rc)
	mF, _, err := mapSampleNames(mapFile)
	if err != nil {
		err = fmt.Errorf("Datachecks 3: %s", err.Error())
		return
	}

	sampleCountMap := len(mF)
	sampleCountVCF := len(vs)

	sampleCountROH := len(rs)
	chrCountROH := len(rc)
	chrNamesVCF := vc

	chrNamesROH := ""
	for chr := range rc {
		if strings.HasPrefix(chr, "chr") {
			if chrNamesROH == "" {
				chrNamesROH = "chr"
			} else if chrNamesROH == "num" {
				chrNamesROH = "mix"
			}
		} else {
			chrNamesROH = "num"
		}
	}

	out := fmt.Sprintf("%s chromosome names in ROH file,%d samples in ROH file,%d chromosomes in ROH file,%s chromosome names in vcf,%d sample count vcf,%d sample count map\n", chrNamesROH, sampleCountROH, chrCountROH, chrNamesVCF, sampleCountVCF, sampleCountMap)
	fmt.Println(out)
	if chrNamesROH != chrNamesVCF {
		mapChrs = true
	}

	// where do samples match? vs is vcf samples, rs is roh samples, mF is forwardmap mR is reverse map
	both, first, second := overlap(vs, rs)
	//fmt.Println(both, first, second)
	if both > 0 && first == 0 && second == 0 {
		// no mapping needed
		fmt.Println("No mapping needed")
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
	chromosomei := strings.TrimPrefix(fi[0], "chr")
	starti := fi[1]

	chromosomej := strings.TrimPrefix(fj[0], "chr")
	startj := fj[1]

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

// return numbers in both, in a not b and in b not a
func overlap(a map[string]bool, b map[string]bool) (both, aOnly, bOnly int) {
	both = 0
	aOnly = 0
	bOnly = 0

	for ka := range a {
		if _, ok := b[ka]; ok {
			both++
		} else {
			aOnly++
		}
	}
	for kb := range b {
		if _, ok := a[kb]; !ok {
			bOnly++
		}
	}
	return
}

// shoudl probably change both to interface .. to generalise
// return the subset of a whose keys are in b
func subset(a map[string]bool, b map[string]string) (c map[string]bool) {
	c = make(map[string]bool)
	for k := range a {
		if _, ok := b[k]; ok {
			c[k] = true
		}
	}
	return

}
