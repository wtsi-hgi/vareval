$$namespaces:
  arv: "http://arvados.org/cwl#"
  cwltool: "http://commonwl.org/cwltool#"

cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement

hints:
  ResourceRequirement:
    ramMin: 4000
    coresMin: 4
    tmpdirMin: 1000
  arv:RuntimeConstraints:
    keep_cache: 1024
    outputDirType: keep_output_dir

inputs:
  - id: executable
    type: File   
  - id: ROH_chr
    type: File
  - id: vcf_file
    type: File[]
  - id: sample_mapping
    type: File
    
steps:

  - id: ROH_calc
    scatter:
      - vcf_file  
    run: roh_go.cwl
    in:
      executable: executable
      ROH_chr: ROH_chr
      vcf_file: vcf_file 
      sample_mapping: sample_mapping       
    out: [output1]


outputs:
  - id: calls
    type: File[]
    outputSource: [ROH_calc/output1]


doc: | 
    produce vcfs of het calls in known ROH regions as part of assessment

