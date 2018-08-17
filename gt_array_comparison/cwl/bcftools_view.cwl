class: CommandLineTool
cwlVersion: v1.0
id: bcftools_extract_sample

requirements:
  DockerRequirement:
    dockerPull: mercury/bcftools-1.6:v2

baseCommand:
  - bcftools
  - view

inputs:
  out_type:
    type:
    - 'null'
    - type: enum
      symbols:
      - 'u'
      - 'o'
      - 'z'
      - 'b'
    inputBinding:
      prefix: -O
  region:
    type: string?
    inputBinding:
      prefix: -r
  sample:
    type: string
    inputBinding:
      prefix: -s
  output_filename:
    type: string
    inputBinding:
      prefix: -o
  input_vcf:
    type: File
    inputBinding:
      position: 1

outputs:
  sample_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
