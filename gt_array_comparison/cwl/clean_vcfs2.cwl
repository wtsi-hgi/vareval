class: CommandLineTool
cwlVersion: v1.0
id: Clean_vcfs2

requirements:
  DockerRequirement:
    dockerPull: mercury/bcftools-1.6:v1

baseCommand:
  - bash
  - -c

arguments:
  - >
    bcftools view -r $(inputs.region) -S $(inputs.sample_list) $(inputs.input_vcf) |
    bcftools norm -c $(inputs.check-ref) -f $(inputs.ref) -m $(inputs.multiallelics) |
    bcftools annotate -O$(inputs.out_type) -x $(inputs.remove-annotations) -o $(inputs.output_filename)

inputs:
  out_type:
    type: enum?
    symbols: ['u','o','z','b']
#    inputBinding:
#      prefix: -O
  region:
    type: str?
 #   inputBinding:
 #     prefix: -r
  sample_list:
    type: File?
#    inputBinding:
#      prefix: -s
  check-ref:
    type: enum?
    symbols: ['e', 'w', 'x', 's']
    default: s
#    inputBinding:
#      prefix: -c
  ref:
    type: File?
#   inputBinding:
#      prefix: -f
  multiallelics:
    type: enum?
    symbols: ['+snps', '+indels', '+both', '+any', '-snps', '-indels', '-both', '-any']
    default: +any
#    inputBinding:
#      prefix: -m
  remove-annotations:
    type: str?
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
  output_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)