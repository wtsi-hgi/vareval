class: CommandLineTool
cwlVersion: v1.0
id: bcftools_extract_sample
baseCommand:
  - bcftools
  - view

inputs:
  out_type:
    type: enum?
    symbols: ['u','o','z','b']
    inputBinding:
      prefix: -O
    doc: Output blah
  sample:
    type: string
    inputBinding:
      prefix: -s
    doc: Blah
  output_filename:	
    type: string
    inputBinding:
      prefix: -o
      valueFrom: "$(runtime.outdir)/$(self)"
  input_vcf:
    type: File
    inputBinding:
      position: 1
    doc: Blah

outputs:
  output_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
