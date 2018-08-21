cwlVersion: v1.0
class: CommandLineTool
doc: Makes an array from the lines contained in a file

requirements:
  - class: InitialWorkDirRequirement
    listing:
    - entry: $(inputs.input_file)
      entryname: $(inputs.input_file.basename)
  - class: DockerRequirement
    dockerPull: python:3

baseCommand:
  - python
  - -c
  - |
      import json
      import sys

      input_filename = sys.argv[1]

      with open(input_filename) as input:
        lines = [l.strip() for l in input.readlines()]

      with open("cwl.output.json", "w") as output:
        output.write(json.dumps({
          "contents_array": lines
          }))

arguments: [$(inputs.input_file.basename)]

inputs:
  - id: input_file
    type: File
    inputBinding:
      loadContents: true

outputs:
  - id: contents_array
    type: string[]
