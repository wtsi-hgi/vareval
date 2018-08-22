cwlVersion: v1.0
class: CommandLineTool
id: Split_vcf_by_af

requirements:
  DockerRequirement:
    dockerPull: mercury/bcftools-1.6:v2

baseCommand:
  - bash

inputs:
  script:
    type: File
    inputBinding:
      position: -1
  input_vcf:
    type: File
    inputBinding:
      position: 1
    doc: Input VCF
  ac_file:
    type: File
    inputBinding:
      position: 2
    doc: AC file, containing allele counts to be split by
  af_file:
    type: File
    inputBinding:
      position: 3
    doc: AF file, containing allele frequencies to be stratified by

outputs:
  stratified_vcfs:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.vcf.gz"
