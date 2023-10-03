cwlVersion: v1.2
class: CommandLineTool
id: phlat_bowtie2
doc: "A specific implementation of bowtie2 to run PHLAT HLA typing"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.threads)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/brownm28/phlat:1.0.0

baseCommand: [tar, xzf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      && bowtie2
  - position: 2
    shellQuote: false
    valueFrom: >-
      | grep -P "\t(A\*|B\*|C\*|DQA1\*|DQB1\*|DRB1\*)" | grep -P -v "\tchr" > $(inputs.output_basename).bwt2.preprocess.sam

inputs:
  orientation: {type: ['null', {type: enum, name: orientation, symbols: ["--fr", "--rf", "--ff"]}], default: "--fr", doc: "If paired-end, read orientation",
    inputBinding: {position: 1}}
  maxins: { type: 'int?', doc: "maximum fragment length", default: 450, inputBinding: { position: 1, prefix: "--maxins" } }
  no_mixed: { type: 'boolean?', doc: "suppress unpaired alignments for paired reads", default: true, inputBinding: { position: 1, prefix: "--no-mixed" } }
  no_discordant: { type: 'boolean?', doc: "suppress discordant alignments for paired reads", default: true, inputBinding: { position: 1, prefix: "--no-discordant" } }
  no_unal: { type: 'boolean?', doc: "suppress unaligned reads", default: true, inputBinding: { position: 1, prefix: "--no-unal" } }
  threads: { type: 'int?', doc: "Num processing threads tp use", default: 8, inputBinding: { position: 1, prefix: "-p" } }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 16 }
  extend_fail: { type: 'int?', doc: "give up extending after <int> failed extends in a row", default: 20, inputBinding: { position: 1, prefix: "-D"} }
  repetitive_seeds: { type: 'int?', doc: "for reads w/ repetitive seeds, try <int> sets of seeds", default: 3, inputBinding: { position: 1, prefix: "-R"} }
  max_seed_mm: { type: 'int?', doc: "max # mismatches in seed alignment; can be 0 or 1", default: 0, inputBinding: { position: 1, prefix: "-N"} }
  seed_substr_len: { type: 'int?', doc: "length of seed substrings; must be >3, <32", default: 20, inputBinding: { position: 1, prefix: "-L"} }
  seed_intervals: { type: 'string?', doc: "interval between seed substrings w/r/t read len (S,1,1.15)", default: "S,1,0.50", inputBinding: { position: 1, prefix: "-i", shellQuote: false} }
  bwt2_index: { type: File, doc: "bowtie2 tar-gzipped index", inputBinding: { position: 0 } }
  bwt2_index_path_prefix: { type: string, doc: "unzipped path + prefix of index, i.e. index4phlat/ucsc.artHLA", inputBinding: { position: 1, prefix: "-x" } }
  unpaired_reads: { type: 'File?', doc: "fastq reads if NOT paired", inputBinding: { position: 1 , prefix: "-U" } }
  read1: { type: 'File?', doc: "read1 fastq reads if paired-end", inputBinding: { position: 1 , prefix: "-1" } }
  mates: { type: 'File?', doc: "mates (read2) fastq reads if paired-end", inputBinding: { position: 1 , prefix: "-2" } }
  output_basename: { type: string, doc: "String to use as basename for output file" }

outputs:
  preprocess_sam: { type: File, outputBinding: { glob: '*.sam' } }
