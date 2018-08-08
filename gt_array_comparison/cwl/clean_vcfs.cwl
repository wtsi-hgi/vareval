cwlVersion: v1.0
class: CommandLineTool
id: Clean_vcfs

requirements:
  DockerRequirement:
    dockerPull: mercury/bcftools-1.6:v2

inputs:
  script:
    type: File
    inputBinding:
      position: -1
    doc:
     Script to preprocess VCF files so they can be fed into hap.py.
     Splits by chrom, sample, normalises REF mismatches, merges multiallelics and removes the AD field.

  ref_file:
    type: File
    inputBinding:
      prefix: -r
    doc:
      Reference genome in fasta format

  sample_list:
    type: File
    inputBinding: 
      prefix: -s
    doc:
      File containing list of samples to be extracted from VCF, one per line

  exome_mask:
    type: File?
    inputBinding:
      prefix: -e
    doc: 
      BED file containing coordinates for exome regions

  input_vcf:
    type: File
    inputBinding:
      position: 1
    doc: Input multi-sample VCF	  
  

outputs:
  vcf_files:
    type: 
      type: array
      item: File
    outputBinding:
      glob: *.vcf.gz
  tbi_files:
    type:
      type: array
      item: File
    outputBinding:
      glob: *.vcf.gz.tbi    
