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
  - id: regions_array
    type: string[]
  - id: multisample_vcf
    type: File
  - id: threads
    type: int?

steps:
  - id: Clean_vcfs2
    run: ./clean_vcfs2.cwl
    scatter: region
    in:
      region: regions_array
      out_type:
        valueFrom: "z"
      ref: fasta_ref
      sample_list: sample_list
      output_filename:
        valueFrom: $(inputs.input_vcf.basename.split(".vcf")[0]).$(inputs.region).clean.vcf.gz
      input_vcf: multisample_vcf
    out:
      [clean_vcf_per_chr]
  - id: Make_samples_array
    run: file_contents_to_array.cwl
    in:
      input_file: sample_list
    out:
      [contents_array]
  - id: Extract_sample
    run: ./bcftools_view.cwl
    scatter: 
      - input_vcf
      - sample
    scatterMethod: flat_crossproduct
    in:
      out_type: 
        valueFrom: "z"
      sample: Make_samples_array/contents_array
      output_filename:
        valueFrom: $(inputs.input_vcf.basename.split(".vcf")[0]).$(inputs.sample).vcf.gz
      input_vcf: Clean_vcfs2/clean_vcf_per_chr
    out:
      [sample_vcf]
  - id: Index_tbi
    run: ../../subrepos/arvados-pipelines/cwl/tools/bcftools/bcftools-index-tbi.cwl
    scatter: 
      - vcf
    in: 
      threads: threads
      vcf: Extract_sample/sample_vcf
      output_filename:
        valueFrom: $(inputs.vcf.basename).tbi
    out:
      [index]
  - id: Combine_sample_vcfs_and_tbis
    run: ../../subrepos/arvados-pipelines/cwl/expression-tools/combine_files.cwl
    scatter:
      - main_file
      - secondary_files
    scatterMethod: dotproduct 
    in:
      main_file: Extract_sample/sample_vcf
      secondary_files: Index_tbi/index
    out:
      [file_with_secondary_files]

outputs:
  - id: out_vcfs
    type: File[]
    outputSource: Combine_sample_vcfs_and_tbis/file_with_secondary_files

