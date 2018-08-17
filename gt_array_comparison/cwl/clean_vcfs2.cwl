class: CommandLineTool
cwlVersion: v1.0
id: Clean_vcfs2

requirements:
  DockerRequirement:
    dockerPull: mercury/bcftools-1.6:v2

baseCommand:
  - bash
  - -c

arguments:
  - >
    bcftools view -r $(inputs.region) -S $(inputs.sample_list.basename) $(inputs.input_vcf.basename) |
    bcftools norm -c $(inputs.checkref) -f $(inputs.ref.basename) -m $(inputs.multiallelics) |
    bcftools annotate -O$(inputs.out_type) -x $(inputs.remove_annotations) -o $(inputs.output_filename)

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
#    inputBinding:
#      prefix: -O
  region:
    type: string?
 #   inputBinding:
 #     prefix: -r
  sample_list:
    type: File?
#    inputBinding:
#      prefix: -s
  checkref:
    type:
    - 'null'
    - type: enum
      symbols:
      - 'e'
      - 'w'
      - 'x'
      - 's'
    default: s
#    inputBinding:
#      prefix: -c
  ref:
    type: File?
#   inputBinding:
#      prefix: -f
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
    default: +any
#    inputBinding:
#      prefix: -m
  remove_annotations:
    type: string?
    default: FORMAT/AD
#   inputBinding:
#     prefix: -x
  output_filename:
    type: string?
    default: output.vcf.gz
#    inputBinding:
#      prefix: -o
  input_vcf:
    type: File
#    inputBinding:
#      position: 1

outputs:
  clean_vcf_per_chr:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)