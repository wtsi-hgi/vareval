cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: StepInputExpressionRequirement

inputs:
  - id: script
    type: File
  - id: fasta_ref
    type: File
  - id: sample_list
    type: File
    inputBinding:
      loadContents: True
  - id: exome_bed
    type: File?
  - id: multisample_vcf
    type: File
  - id: threads
    type: int?

steps:
  - id: Clean_vcfs
    run: clean_vcfs.cwl
    in:
      script: script
      ref_file: fasta_ref
      sample_list: sample_list
      exome_mask: exome_bed
      input_vcf: multisample_vcf
    out:
      - chr_vcf_files
      - chr_tbi_files
  - id: Combine_vcfs_and_tbis
    run: ../expression-tools/combine_files.cwl
    scatter:
      - main_file
      - secondary_file
    scatterMethod: dotproduct
    in:
      main_file: Clean_vcfs/chr_vcf_files
      secondary_file: Clean_vcfs/chr_tbi_files
    out: [vcf_with_tbi_files]
  - id: Extract_sample
    run: extract_single_sample_vcf.cwl
    scatter: 
      - sample
      - input_vcf
    scatterMethod: nested_crossproduct
    in:
      out_type: 
        valueFrom: "z"
      sample: 
	valueFrom: $(inputs.sample_list.contents.split("\n"))
      output_filename: $(inputs.input_vcf.basename.split(".vcf")[0]).$(sample).vcf.gz
      input_vcf: Combine_vcfs_and_tbis/vcf_with_tbi_files     
    out:
      - [[sample_chr_vcfs]]
  - id: Flatten-sample-chr-outputs
    run: ../expression-tools/flatten-array-file.cwl
    in:
      2d-array: Extract_sample/sample_chr_vcfs
    out: [flattened_array_sample_chr_vcfs]
  - id: Index_tbi
    run: bcftools_index_tbi.cwl
    scatter: 
      - vcf
    in: 
      threads: threads
      vcf: Flatten-sample-chr-outputs/flattened_array_sample_chr_vcfs
      output_filename: $(inputs.vcf.basename).tbi
    out:
      - [index_files]
  - id: Combine_sample_vcfs_and_tbis
    run: ../expression-tools/combine_files.cwl
    scatter:
      - main_file
      - secondary_file
    scatterMethod: dotproduct 
    in:
      main_file:  Flatten-sample-chr-outputs/flattened_array_sample_chr_vcfs
      secondary_file: Index_tbi/index_files
    out: [sample_vcf_with_tbi_files]

outputs:
  - id: out
    type: File[]
    outputSource: Combine_sample_vcfs_and_tbis/sample_vcf_with_tbi_files