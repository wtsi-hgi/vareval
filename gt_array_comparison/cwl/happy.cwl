class: CommandLineTool
cwlVersion: v1.0
id: run_hap.py

requirements:
  - class: DockerRequirement
    dockerPull: pkrusche/hap.py
  - class: InlineJavascriptRequirement
    expressionLib:
    - |-
      function generate_stratif_file_list(input) {
          if (!input) {
              return "";
          }

          if(!Array.isArray(input)){
              input = [input];
          }

          var output = [];
          input.forEach(function(element) {
              if (element.basename.includes("AF"))
                  var region_name = "AF" + element.basename.split("AF")[1].split(".bed")[0];
              else if (element.basename.includes("AC"))
                  var region_name = "AC" + element.basename.split("AC")[1].split(".bed")[0];
              output.push(region_name + "\t" + element.path);
        })

        return output.join("\n");
      }
  - class: InitialWorkDirRequirement
    listing:
      - entryname: happy.stratif_file
        entry: $(generate_stratif_file_list(inputs.stratification))

baseCommand:
  - /opt/hap.py/bin/hap.py

arguments:
  - "--stratification"
  - "happy.stratif_file"

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
    type:
    - 'null'
    - type: array
      items: File
      inputBinding:
        valueFrom: $(null)
  stratification-fixchr:
    type: boolean?
    inputBinding:
      prefix: --stratification-fixchr
  fixchr:
    type: boolean?
    inputBinding:
      prefix: --fixchr
  sdf_ref_dir:
    type: Directory?
    inputBinding:
      prefix: --engine-vcfeval-template

outputs:
  happy_out:
    type: File[]
    outputBinding:
      glob: $(inputs.output_fileprefix + "*")

