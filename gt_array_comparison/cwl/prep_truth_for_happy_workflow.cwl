cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:
  - id: freq_split_script
    type: File
  - id: truth_vcf
    type: File
  - id: ac_file
    type: File
  - id: af_file
    type: File
  - id: regions
    type: File
  - id: sample_list
    type: File
  - id: fasta_ref
    type: File
  - id: threads
    type: int?

steps:
  - id: Make_regions_array
    run: file_contents_to_array.cwl
    in:
      input_file: regions
    out:
      [contents_array]
  - id: Norm_vcf
    run: bcftools_norm.cwl
    in:
      input_vcf: truth_vcf
      ref: fasta_ref
      out_type:
        valueFrom: "z"
      output_filename:
        valueFrom: $(inputs.input_vcf.basename.split(".vcf")[0]).norm.vcf.gz
      checkref:
        valueFrom: "s"
      multiallelics:
        valueFrom: "+any"
    out:
      - normed_vcf
  - id: Index_multisample_tbi
    run: ../../subrepos/arvados-pipelines/cwl/tools/bcftools/bcftools-index-tbi.cwl
    in:
      threads: threads
      vcf: Norm_vcf/normed_vcf
      output_filename:
        valueFrom: $(inputs.vcf.basename).tbi
    out:
      [index]
  - id: Combine_multisample_vcf_and_tbi
    run: ../../subrepos/arvados-pipelines/cwl/expression-tools/combine_files.cwl
    in:
      main_file:  Norm_vcf/normed_vcf
      secondary_files: Index_multisample_tbi/index
    out:
      [file_with_secondary_files]
  - id: Split_by_af
    run: split_vcf_by_af.cwl
    scatter: region
    in:
      script: freq_split_script
      input_vcf: Norm_vcf/normed_vcf
      ac_file: ac_file
      af_file: af_file
      region: Make_regions_array/contents_array
    out:
      [stratified_beds]
  - id: Make_samples_array
    run: file_contents_to_array.cwl
    in:
      input_file: sample_list
    out:
      [contents_array]
  - id: Extract_region_sample
    run: ./bcftools_view.cwl
    scatter:
      - region
      - sample
    scatterMethod: flat_crossproduct
    in:
      out_type:
        valueFrom: "z"
      sample: Make_samples_array/contents_array
      output_filename:
        valueFrom: $(inputs.input_vcf.basename.split(".vcf")[0]).$(inputs.region).$(inputs.sample).vcf.gz
      input_vcf: Combine_multisample_vcf_and_tbi/file_with_secondary_files
      region: Make_regions_array/contents_array
    out:
      [sample_vcf]
  - id: Index_tbi
    run: ../../subrepos/arvados-pipelines/cwl/tools/bcftools/bcftools-index-tbi.cwl
    scatter:
      - vcf
    in:
      threads: threads
      vcf: Extract_region_sample/sample_vcf
      output_filename:
        valueFrom: $(inputs.vcf.basename).tbi
    out:
      [index]
  - id: Make_conf_regions_beds
    run: extract_regions_from_vcf.cwl
    scatter: vcf
    in:
      vcf: Extract_region_sample/sample_vcf
      output_filename:
        valueFrom: $(inputs.vcf.basename.split(".vcf")[0]).bed
    out:
      [bed]
  - id: Combine_sample_vcfs_and_tbis
    run: ../../subrepos/arvados-pipelines/cwl/expression-tools/combine_files.cwl
    scatter:
      - main_file
      - secondary_files
    scatterMethod: dotproduct
    in:
      main_file: Extract_region_sample/sample_vcf
      secondary_files: Index_tbi/index
    out:
      [file_with_secondary_files]


outputs:
  - id: out_vcfs
    type: File[]
    outputSource: Combine_sample_vcfs_and_tbis/file_with_secondary_files
  - id: truth_beds
    type: File[]
    outputSource: Make_conf_regions_beds/bed
  - id: af_stratif_beds
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: Split_by_af/stratified_beds
