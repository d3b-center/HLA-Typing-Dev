cwlVersion: v1.2
class: CommandLineTool
id: t1k_genotyper
doc: "Run T1K genotyper 'The ONE genotyper for Kir and HLA'"
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: $(inputs.ram * 1000)
    coresMin: $(inputs.threads)
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/brownm28/t1k:v1.0.2
baseCommand: [tar xzf]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      && perl /T1K-1.0.2/run-t1k

inputs:
  reads: { type: 'File?', doc: "read1 fastq reads if paired-end", inputBinding: { position: 1 , prefix: "-1" } }
  mates: { type: 'File?', doc: "mates (read2) fastq reads if paired-end", inputBinding: { position: 2 , prefix: "-2" } }
  single_end: { type: 'File?', doc: "fastq reads if single-end", inputBinding: { position: 1 , prefix: "-u" } }
  interleaved: { type: 'File?', doc: "fastq reads if interleaved", inputBinding: { position: 1 , prefix: "-i" } }
  bam: { type: 'File?', doc: "if bam input", inputBinding: { position: 1 , prefix: "-n" } }
  t1k_ref_tar: { type: File, doc: "Tar ball with T1K reference", inputBinding: { position: 0 }}
  t1k_ref_path: { type: string, doc: "String of path to un-tarred hla ref, i.e. kiridx/kiridx_rna_seq.fa", inputBinding: { position: 2, prefix: "-f"}  }
  preset: {type: ['null', {type: enum, name: preset, symbols: ["hla", "hla-wgs", "kir-wgs", "kir-wes"]}], default: "hla", doc: "If paired-end, read orientation",
    inputBinding: { position: 2 , prefix: "--preset"} }
  threads: { type: 'int?', doc: "Num processing threads to use", default: 8, inputBinding: { position: 2, prefix: "-t" } }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 16 }
  output_basename: { type: 'string?', doc: "Prefix string for output file names. Default inferred from input", inputBinding: { position: 2, prefix: "-o"} }
outputs:
  aligned_fasta: { type: 'File[]', outputBinding: { glob: '*_aligned_*.fa' } }
  allele_tsv: { type: File, outputBinding: { glob: '*_allele.tsv' } }
  allele_vcf: { type: File, outputBinding: { glob: '*_allele.vcf' } }
  candidate_fastqs: { type: 'File[]', outputBinding: { glob: '*_candidate_*.fq' } }
  genotype_tsv: { type: File, outputBinding: { glob: '*_genotype.tsv' } }