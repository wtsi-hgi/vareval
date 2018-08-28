class: CommandLineTool
cwlVersion: v1.0
id: run_hap.py

requirements:
  DockerRequirement:
    dockerPull: pkrusche/hap.py

baseCommand:
  - /opt/hap.py/libexec/rtg-tools-install/rtg
  - format

inputs:
  ref_fasta:
    type: File
    inputBinding:
      position: 1
  output_dirname:
    type: string?
    inputBinding:
      prefix: -o

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dirname)