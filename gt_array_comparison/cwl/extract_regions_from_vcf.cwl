cwlVersion: v1.0
class: CommandLineTool
id: Make_bed_from_vcf

requirements:
  DockerRequirement:
    dockerPull: ubuntu:14.04
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.vcf)
        entryname: $(inputs.vcf.basename)

baseCommand:
  - "bash"
  - "-c"

arguments:
  - zcat $(inputs.vcf.basename) | grep -v "^#" | awk 'BEGIN { OFS="\\t"; } { print $1, $2-1, $2-1+length($4); }' > $(inputs.output_filename)

inputs:
  vcf:
    type: File
  output_filename:
    type: string?
    default: output.bed

outputs:
  bed:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)