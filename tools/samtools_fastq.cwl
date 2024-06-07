cwlVersion: v1.2
class: CommandLineTool
id: samtools_fastq
doc: |
  Convert SAM to fastq
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'pgc-images.sbgenomics.com/d3b-bixu/samtools:1.20-parallel'
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 32000
    coresMin: $(inputs.threads)

baseCommand: []
arguments:
  - position: 10
    shellQuote: false
    valueFrom: >-
      samtools collate -@ $(inputs.threads - 1) -u -O $(inputs.input_reads.path) | samtools fastq -s /dev/null -0 /dev/null -n > $(inputs.output_basename != null ? inputs.output_basename + ".fastq" : inputs.input_reads.basename.replace(/.bam$/,".fastq"))

inputs:
  input_reads: {type: File, secondaryFiles: ['.bai']}
  threads:
    type: int? 
    default: 16
  output_basename: {type: 'string?', doc: "String to use as basename for outputs"}
outputs:
  fastq:
    type: File
    outputBinding:
      glob: '*.fastq'
