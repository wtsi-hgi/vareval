cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: ScatterFeatureRequirement
  - class: InlineJavascriptRequirement


inputs:
  - id: fasta_ref
    type: File
  - id: sample_list
    type: File
  - id: regions
    type: File?
  - id: test_vcf
    type: File
  - id: truth_vcf
    type: File
  - id: ac_file
    type: File
  - id: af_file
    type: File
  - id: freq_split_script
    type: File
  - id: index_threads
    type: int?

steps:
  - id: Prep_truth
    run: prep_truth_for_happy_workflow.cwl
    in:
      freq_split_script: freq_split_script
      fasta_ref: fasta_ref
      truth_vcf: truth_vcf
      ac_file: ac_file
      af_file: af_file
      sample_list: sample_list
      regions: regions
    out:
      - out_vcfs
      - truth_beds
      - af_stratif_beds
  - id: Prep_test
    run: prep_vcfs_for_happy_workflow.cwl
    in:
      fasta_ref: fasta_ref
      sample_list: sample_list
      regions: regions
      multisample_vcf: test_vcf
    out:
      [out_vcfs]
  - id: Prep_ref
    run: format_sdf_ref.cwl
    in:
      ref_fasta: fasta_ref
      output_dirname:
        valueFrom: $(inputs.ref_fasta.basename.split(".fa")[0]).sdf
    out:
      - output_dir
  - id: Format_stratification_beds
    run: multiply_inner_arrays_in_2d_arrays.cwl
    in:
      2d_array: Prep_truth/af_stratif_beds
      multiplication_factor:
        valueFrom: $( 3 )
    out:
      - multiplied-2d-array
  - id: Run_happy
    run: happy.cwl
    scatter:
      - truth_vcf
      - test_vcf
      - confident_regions
      - stratification_beds
    scatterMethod: dotproduct
    in:
      truth_vcf: Prep_truth/out_vcfs
      test_vcf: Prep_test/out_vcfs
      confident_regions: Prep_truth/truth_beds
      ref: fasta_ref
      sdf_ref_dir: Prep_ref/output_dir
      stratification_beds: Format_stratification_beds/multiplied-2d-array
      output_fileprefix:
        valueFrom: $("output." + inputs.truth_vcf.basename.split(".vcf")[0].split(".").slice(-2).join("."))
    out:
      [happy_out]

outputs:
  - id: eval_out
    type:
      type: array
      items:
        type: array
        items: File
    outputSource: Run_happy/happy_out


