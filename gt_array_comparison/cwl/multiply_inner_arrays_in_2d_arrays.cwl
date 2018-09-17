cwlVersion: v1.0
class: CommandLineTool

requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: inputs.json
        entry: '{"2d_array": $(inputs.2d_array)}'
  DockerRequirement:
    dockerPull: python:3

baseCommand:
  - python
  - -c
  - |
      import json
      import sys

      factor = int(sys.argv[1])

      with open("inputs.json") as inputs_file:
        input_matrix = json.load(inputs_file)["2d_array"]

      output_matrix= []

      for inner_array in input_matrix:
        output_matrix.extend([inner_array] * factor)

      with open("cwl.output.json", "w") as output_file:
        output_file.write(json.dumps({"multiplied-2d-array": output_matrix}))

arguments: [$(inputs.multiplication_factor)]

inputs:
  2d_array:
    type:
      type: array
      items:
        type: array
        items: File
  multiplication_factor:
    type: int

outputs:
  multiplied-2d-array:
    type:
      type: array
      items:
        type: array
        items: File