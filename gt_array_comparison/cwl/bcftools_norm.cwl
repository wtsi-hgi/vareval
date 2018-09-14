class: CommandLineTool
cwlVersion: v1.0
id: bcftools_norm

requirements:
  DockerRequirement:
    dockerPull: mercury/bcftools-1.6:v2

baseCommand:
  - bcftools
  - norm

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
  checkref:
    type:
    - 'null'
    - type: enum
      symbols:
      - 'e'
      - 'w'
      - 'x'
      - 's'
    inputBinding:
      prefix: -c
  ref:
    type: File?
    inputBinding:
      prefix: -f
  multiallelics:
    type:
    - 'null'
    - type: enum
      symbols:
      - '+snps'
      - '+indels'
      - '+both'
      - '+any'
      - '-snps'
      - '-indels'
      - '-both'
      - '-any'
    inputBinding:
      prefix: -m
  output_filename:
    type: string?
    default: output.vcf.gz
    inputBinding:
      prefix: -o
  input_vcf:
    type: File
    inputBinding:
      position: 1

outputs:
  normed_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)