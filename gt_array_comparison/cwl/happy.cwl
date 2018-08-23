class: CommandLineTool
cwlVersion: v1.0
id: run_hap.py

requirements:
  DockerRequirement:
    dockerPull: pkrusche/hap.py

baseCommand:
  - /opt/hap.py/bin/hap.py

inputs:
  truth_vcf:
    type: File
    inputBinding:
      position: 1
  test_vcf:
    type: File
    inputBinding:
      position: 2
  ref:
    type: File
    inputBinding:
      prefix: -r
  write_vcf:
    type: boolean?
    inputBinding:
      prefix: --write-vcf
    default: true
  output_fileprefix:
    type: string?
    inputBinding:
      prefix: -o
    default: output
  eval_engine:
    type:
    - 'null'
    - type: enum
      symbols:
      - 'vcfeval'
      - 'xcmp'
      - 'scmp-distance'
    inputBinding:
      prefix: --engine
    default: 'vcfeval'
  confident_regions:
    type: File
    inputBinding:
      prefix: -f
  filter:
    type: boolean?
    inputBinding:
      prefix: --pass-only
    default: true
  stratification:
    type: File?
    inputBinding:
      prefix: --stratification
  sdf_ref_name:
    type: string?
    inputBinding:
      prefix: --engine-vcfeval-template

outputs:
  happy_stdout:
    type: stdout
  happy_out:
    type: File[]
    outputBinding:
      glob: $(inputs.output_fileprefix)

