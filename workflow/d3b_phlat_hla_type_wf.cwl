cwlVersion: v1.2
class: Workflow
id: d3b-phlat-hla-type-wf
doc: "D3b implementation of PHLAT HLA Typing: https://pubmed.ncbi.nlm.nih.gov/29858810/"
requirements:
- class: SubworkflowFeatureRequirement
inputs:
  # common
  output_basename: { type: string, doc: "String to use as basename for output file" }
  # Bowtie2
  orientation: {type: ['null', {type: enum, name: orientation, symbols: ["--fr", "--rf", "--ff"]}], default: "--fr", doc: "If paired-end, read orientation" }
  maxins: { type: 'int?', doc: "maximum fragment length", default: 450 }
  no_mixed: { type: 'boolean?', doc: "suppress unpaired alignments for paired reads", default: true }
  no_discordant: { type: 'boolean?', doc: "suppress discordant alignments for paired reads", default: true }
  no_unal: { type: 'boolean?', doc: "suppress unaligned reads", default: true }
  bwt2_threads: { type: 'int?', doc: "Num processing threads to use for bowtie2", default: 8 }
  bwt2_ram: { type: 'int?', doc: "Num GB memory to make available for bowtie2", default: 16 }
  extend_fail: { type: 'int?', doc: "give up extending after <int> failed extends in a row", default: 20 }
  repetitive_seeds: { type: 'int?', doc: "for reads w/ repetitive seeds, try <int> sets of seeds", default: 3 }
  max_seed_mm: { type: 'int?', doc: "max # mismatches in seed alignment; can be 0 or 1", default: 0 }
  seed_substr_len: { type: 'int?', doc: "length of seed substrings; must be >3, <32", default: 20 }
  seed_intervals: { type: 'string?', doc: "interval between seed substrings w/r/t read len (S,1,1.15)", default: "S,1,0.50" }
  bwt2_index: { type: File, doc: "bowtie2 tar-gzipped index" }
  bwt2_index_path_prefix: { type: string, doc: "unzipped path + prefix of index, i.e. index4phlat/ucsc.artHLA" }
  unpaired_reads: { type: 'File?', doc: "fastq reads if NOT paired" }
  reads: { type: 'File?', doc: "read1 fastq reads if paired-end" }
  mates: { type: 'File?', doc: "mates (read2) fastq reads if paired-end" }
  # Format sam
  bwt2_header: { type: File, doc: "Desired header file to prepend to processed sorted bam" }
  # phlat scoring
  preload_pickle: { type: File, doc: "pickle file from PHLAT resources" }
  
outputs:
  phlat_hla_summary: {type: File, outputSource: phlat_hla_scoring/phlat_HLA_summary}
  phlat_qnsorted_bam: {type: File, outputSource: phlat_format_sam/phlat_name_sort}
  phlat_hla_map: {type: File, outputSource: phlat_format_sam/hla_map}
steps:
  phlat_bowtie2:
    run: ../tools/phlat_bowtie2.cwl
    in:
      orientation: orientation
      maxins: maxins
      no_mixed: no_mixed
      no_discordant: no_discordant
      no_unal: no_unal
      threads: bwt2_threads
      ram: bwt2_ram
      extend_fail: extend_fail
      repetitive_seeds: repetitive_seeds
      max_seed_mm: max_seed_mm
      seed_substr_len: seed_substr_len
      seed_intervals: seed_intervals
      bwt2_index: bwt2_index
      bwt2_index_path_prefix: bwt2_index_path_prefix
      unpaired_reads: unpaired_reads
      reads: reads
      mates: mates
      output_basename: output_basename
    out: [preprocess_sam]
  phlat_format_sam:
    run: ../tools/phlat_format_sam.cwl
    in:
      phlat_sam: phlat_bowtie2/preprocess_sam
      bwt2_header: bwt2_header
    out: [hla_map, phlat_bam, phlat_name_sort]
  phlat_hla_scoring:
    run: ../tools/phlat_hla_scoring.cwl
    in:
      qsorted_bam: phlat_format_sam/phlat_name_sort
      map_file: phlat_format_sam/hla_map
      preload_pickle: preload_pickle
      output_basename: output_basename
    out: [phlat_HLA_summary]

$namespaces:
  sbg: https://sevenbridges.com

