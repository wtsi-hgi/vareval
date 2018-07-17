#$$namespaces:
#  arv: "http://arvados.org/cwl#"
#  cwltool: "http://commonwl.org/cwltool#"

cwlVersion: v1.0
class: Workflow

requirements:
  - class: MultipleInputFeatureRequirement

hints:
  ResourceRequirement:
    ramMin: 4000
    coresMin: 4
    tmpdirMin: 1000
#  arv:RuntimeConstraints:
##    keep_cache: 1024
#    outputDirType: keep_output_dir

inputs:
  - id: executable
    type: File   
  - id: ROH_chr
    type: File
  - id: vcf_file
    type: File
  - id: sample_mapping
    type: File
    
steps:

# check whether sample and chromosome namings match in vcf and roh file
# 
- id: data_checks
    run: data_checks.cwl
    in:
      executable: executable
      ROH_chr: ROH_chr
      vcf_file: vcf_file 
      sample_mapping: sample_mapping       
    out: [output1, output2]  # map chromosomes, map samples

# make individual bed files per sample from roh file
  - id: make_sample_beds
    run: sample_beds.cwl
    in:
      executable: executable
      ROH_chr: ROH_chr
      vcf_file: vcf_file 
      sample_mapping: sample_mapping       
    out: [output1, output2]  # bed files and log file

# 

outputs:
  - id: calls
    type:
      type: array
      items:
        - type: array
          items: File
    outputSource: [ROH_calc/output1]
  - id: stats
    type: File[]
    outputSource: [ROH_calc/output2]


doc: | 
    produce vcfs of het calls in known ROH regions as part of vareval assessment

