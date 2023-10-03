cwlVersion: v1.2
class: CommandLineTool
id: phlat_hla_scoring
doc: "D3b rewrite of PHLAT.py wrapper script"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.threads)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/brownm28/phlat:1.0.0
baseCommand: [python]
arguments:
  - position: 0
    shellQuote: false
    valueFrom: >-
      /HLA-Typing-Dev/scripts/phlat_d3b.py

inputs:
  qsorted_bam: { type: File, doc: "Queryname sorted bam", inputBinding: { position: 0, prefix: "--qsorted-bam" } }
  map_file: { type: File, doc: "HLA map file", inputBinding: { position: 0, prefix: "--map-file" } }
  preload_pickle: { type: string, doc: "pickle file from PHLAT resources", inputBinding: { position: 0, prefix: "--preload-pickle" } }
  output_basename: { type: string, doc: "String to use as basename for output file", inputBinding: { position: 0, prefix: "--output-basename" } }
  threads: { type: 'int?', default: 4 }
  ram: { type: 'int?', default: 4 }
outputs:
  phlat_HLA_summary: { type: File, outputBinding: { glob: '*_HLA.sum' } }
