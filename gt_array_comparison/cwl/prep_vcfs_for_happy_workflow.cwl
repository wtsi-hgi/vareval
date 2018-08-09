cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:
  - id: fasta_ref
    type: File
  - id: sample_list
    type: File
    inputBinding:
      loadContents: True
  - id: regions
    type: File
    inputBinding:
      loadContents: True
  - id: multisample_vcf
    type: File
  - id: threads
    type: int?

steps:
  - id: Clean_vcfs2
    run: ./clean_vcfs2.cwl
    scatter: region
    in:
      region:
        source: regions
        valueFrom: $(self.contents.split("\n"))
      out_type:
        valueFrom: "z"
      ref: fasta_ref
      sample_list: sample_list
      output_filename: $(inputs.input_vcf.basename.split(".vcf")[0]).$(inputs.region).clean.vcf.gz
      input_vcf: multisample_vcf
    out: chr_vcf_files
  - id: Extract_sample
    run: ./extract_single_sample_vcf.cwl
    scatter: 
      - sample
      - input_vcf
    scatterMethod: flat_crossproduct
    in:
      out_type: 
        valueFrom: "z"
      sample:
        source: sample_list
    	valueFrom: $(self.contents.split("\n"))
      output_filename: $(inputs.input_vcf.basename.split(".vcf")[0]).$(sample).vcf.gz
      input_vcf: Clean_vcfs2/chr_vcf_files
    out:
      - sample_chr_vcfs
  - id: Index_tbi
    run: ../subrepos/arvados-pipelines/cwl/bcftools_index_tbi.cwl
    scatter: 
      - vcf
    in: 
      threads: threads
      vcf: Extract_sample/sample_chr_vcfs
      output_filename: $(inputs.vcf.basename).tbi
    out:
      - index_files
  - id: Combine_sample_vcfs_and_tbis
    run: ../subrepos/arvados-pipelines/cwl/expression-tools/combine_files.cwl
    scatter:
      - main_file
      - secondary_file
    scatterMethod: dotproduct 
    in:
      main_file: Extract_sample/sample_chr_vcfs
      secondary_file: Index_tbi/index_files
    out: sample_vcf_with_tbi_files

outputs:
  - id: out
    type: File[]
    outputSource: Combine_sample_vcfs_and_tbis/sample_vcf_with_tbi_files

