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
  preset: {type: ['null', {type: enum, name: preset, symbols: ["hla", "hla-wgs", "kir-wgs", "kir-wes"]}], default: "hla", doc: "preset parameters for cases requiring non-default settings",
    inputBinding: { position: 2 , prefix: "--preset"} }
  threads: { type: 'int?', doc: "Num processing threads to use", default: 8, inputBinding: { position: 2, prefix: "-t" } }
  ram: { type: 'int?', doc: "Num GB memory to make available", default: 16 }
  output_basename: { type: 'string?', doc: "Prefix string for output file names. Default inferred from input", inputBinding: { position: 2, prefix: "-o"} }
  gene_coords: { type: 'File?', doc: "path to the gene coordinate file (required when -b input)",
    inputBinding: { position: 1, prefix: "-c" } }
  min_align_similarity: { type: 'float?', doc: "minimum alignment similarity", default: 0.8,
    inputBinding: { position: 1, prefix: "-s" } }
  frac: { type: 'float?', doc: "filter if abundance is less than the frac of dominant allele", default: 0.15,
    inputBinding: { position: 1, prefix: "--frac" } }
  cov: { type: 'float?', doc: "filter genes with average coverage less than the specified value", default: 1.0,
    inputBinding: { position: 1, prefix: "--cov"} }
  crossGeneRate: { type: 'float?', doc: "the effect from other gene's expression", default: 0.04,
    inputBinding: { position: 1, prefix: "--crossGeneRate"} }
  alleleDigitUnits: { type: 'int?', doc: "the number of units in genotyping result. default value is auto-determined by algorthim",
    inputBinding: { position: 1, prefix: "--alleleDigitUnits" } }
  alleleDelimiter: { type: 'string?',  doc: "the delimiter character for digit unit.default value is auto-determined by algorthim",
    inputBinding: { position: 1, prefix: "--alleleDelimiter" } }
  alleleWhitelist: { type: 'string?', doc: "only consider read aligned to the listed allele series. not used by default",
    inputBinding: { position: 1 ,prefix: "--alleleWhitelist" } }
  barcode: { type: 'File?', doc: "if -b, BAM field for barcode; if -1 -2/-u, file containing barcodes",
    inputBinding: { position: 1, prefix: "--barcode"} }
  barcodeRange: { type: 'string[]?', doc: "INT INT CHAR: start, end(-1 for length-1), strand in a barcode is the true barcode (default: 0 -1 +)",
    inputBinding: { position: 1, prefix: "--barcodeRange" } }
  barcodeWhitelist: { type: 'File?', doc: "path to the barcode whitelist (default: not used)",
    inputBinding: { position: 1, prefix: "--barcodeWhitelist" } }
  read1Range: { type: 'int[]?', doc: "start, end(-1 for length-1) in -1/-u files for genomic sequence (default: 0 -1)",
    inputBinding: { position: 1, prefix: "--read1Range" } }
  read2Range: { type: 'int[]?', doc: "start, end(-1 for length-1) in -2 files for genomic sequence (default: 0 -1)",
    inputBinding: { position: 1, prefix: "--read2Range" } }
  mateIdSuffixLen: { type: 'int?', doc: "the suffix length in read id for mate. (default: not used)",
    inputBinding: { position: 1, prefix: "--mateIdSuffixLen" } }
  abnormalUnmapFlag: { type: 'string?', doc: "the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)",
    inputBinding: { position: 1, prefix: "--abnormalUnmapFlag" } }
  relaxIntronAlign: { type: 'boolean?', doc: "allow one more mismatch in intronic alignment",  default: false,
    inputBinding: { position: 1, prefix: "--relaxIntronAlign" } }
  noExtraction: { type: 'boolean?', doc: "directly use the files from provided -1 -2/-u for genotyping", default: false,
    inputBinding: { position: 1, prefix: "--noExtraction"} }
  skipPostAnalysis: { type: 'boolean?', doc: "only conduct genotyping.", default: false,
    inputBinding: { position: 1, prefix: "--skipPostAnalysis"} }
outputs:
  aligned_fasta: { type: 'File[]', outputBinding: { glob: '*_aligned_*.fa' } }
  allele_tsv: { type: File, outputBinding: { glob: '*_allele.tsv' } }
  allele_vcf: { type: File, outputBinding: { glob: '*_allele.vcf' } }
  candidate_fastqs: { type: 'File[]', outputBinding: { glob: '*_candidate_*.fq' } }
  genotype_tsv: { type: File, outputBinding: { glob: '*_genotype.tsv' } }