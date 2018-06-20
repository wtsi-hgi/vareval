The aim is to find heterozygous calls that are in known homozygous regions.

This is part of evaluating different variant callers, for the ELGH dataset
which has the known ROH data.

There are two complications:

1) the sample names differ in the multisample vcfs 
and ROH data, so one input is a mapping file between the two sample names. 

2) The chromosome names are 1 in the ROH data and chr1 in the vcfs. 

Inputs are 

1) an array of multisample vcfs that we have ROH regions for
2) the sample id mapping file
3) the ROH regions per chromosome and sample (one file)

Outputs are
vcfs of any het calls in a known ROH region, one file per input vcf and sample
named: het_in_roh_<inputvcf>_<sampleid>.vcf.gz

The program is in Go because I wanted to experiment with loading an
executable in cwl and also becasue it's nicer for error handling than 
bash but basically I'm using it as a script to link bcftools calls 
and it could easily be changed to python or bash

It has run OK in Arvados with all inputs in keep and a yaml file like this:
--------------------------------------------------------------------------------
ROH_chr:
  class: File
  location:  keep:<content address>/allROH.txt
   
executable:
  class: File
  location: keep:<content address>/roh_comparison

vcf_file:
  - class: File
    location: keep:<content address>/ELGH-V2plus.gatk_gvcf_to_vcf.20180516.hard-filtered.vcf.gz
    secondaryFiles:
      - class: File
        location:  keep:<content address>/ELGH-V2plus.gatk_gvcf_to_vcf.20180516.hard-filtered.vcf.gz.tbi
  - class: File
    location: keep:<content address>/ELGH-V2plus.gatk_gvcf_to_vcf.20180516.vcf.gz
    secondaryFiles:
      - class: File
        location:  keep:<content address>/ELGH-V2plus.gatk_gvcf_to_vcf.20180516.vcf.gz.tbi
  - class: File
    location: keep:<content address>/ELGH-V2plus.gvcf_to_vcf.20180516.hard-filtered.vcf.gz
    secondaryFiles:
      - class: File
        location:  keep:<content address>/ELGH-V2plus.gvcf_to_vcf.20180516.hard-filtered.vcf.gz.tbi
  - class: File
    location: keep:<content address>/ELGH-V2plus.gvcf_to_vcf.20180516.vcf.gz
    secondaryFiles:
      - class: File
        location:  keep:<content address>/ELGH-V2plus.gvcf_to_vcf.20180516.vcf.gz.tbi

sample_mapping:
  class: File
  location: keep:<content address>/sample_id_mappings_egan_to_elgh
----------------------------------------------------------------------------- 
  and command line (something like)
  
   arvados-cwl-runner --api=containers --no-log-timestamps --disable-reuse 
   --on-error=continue --eval-timeout=20 --output-name=sarah --output-tags=tag1
    roh_go_workflow.cwl roh_go_workflow.yaml 

