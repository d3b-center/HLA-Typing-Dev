cwlVersion: v1.2
class: Workflow
id: d3b-phlat-hla-type-wf
doc: |+
  # D3b PHLAT Implementation

  This is a software revival of sorts based on the the work from [PHLAT: Inference of High-Resolution HLA Types from RNA and Whole Exome Sequencing](https://pubmed.ncbi.nlm.nih.gov/29858810/).
  It was unfortunately implemented in python 2.7 and currently there are no plans to attempt to port over the original code to python3.
  However, the workflow implementation should make this far easier and more accessible to use.
  Some of the support scripts, like `scripts/extensions.py` and `scripts/utilities.py` have been copied over as-is to help ensure functional equivalency, while the former `PHLAT.py` (`go()` function replaced by `scripts/phlat_d3b.py`) and `prepare.py` have been replaced with cwltool implementations for ease of use and readability

  ## [D3b PHLAT Workflow](workflow/d3b_phlat_hla_type_wf.cwl)

  ### Inputs
  Common:
   - `output_basename`: String to use as basename for output file
  Required for bowtie2:
   - `bwt2_index`: bowtie2 tar-gzipped index
   - `bwt2_index_path_prefix`: unzipped path + prefix of index, i.e. index4phlat/ucsc.artHLA
   - `unpaired_reads`: fastq reads if NOT paired
   - `reads`: read1 fastq reads if paired-end
   - `mates`: mates (read2) fastq reads if paired-end
  Required for SAM output formatting:
   - `bwt2_header`: Desired header file to prepend to processed sorted bam
  Required for final PHLAT HLA Scoring
   - `preload_pickle`: pickle file from PHLAT resources
  Optional bowtie2:
   - `orientation`: If paired-end, read orientation. Default value: "--fr"
   - `maxins`: maximum fragment length. Default value: 450
   - `no_mixed`: suppress unpaired alignments for paired reads. Default value: true
   - `no_discordant`: suppress discordant alignments for paired reads. Default value: true
   - `no_unal`: suppress unaligned reads. Default value: true
   - `bwt2_threads`: Num processing threads to use for bowtie2. Default value: 8
   - `bwt2_ram`: Num GB memory to make available for bowtie2. Default value: 16
   - `extend_fail`: give up extending after <int> failed extends in a row. Default value: 20
   - `repetitive_seeds`: for reads w/ repetitive seeds, try <int> sets of seeds. Default value: 3
   - `max_seed_mm`: max # mismatches in seed alignment; can be 0 or 1. Default value: 0
   - `seed_substr_len`: length of seed substrings; must be >3, <32. Default value: 20
   - `seed_intervals`: interval between seed substrings w/r/t read len (S,1,1.15). Default value: "S,1,0.50"
  ### Outputs
   - `phlat_hla_summary`: Final PHLAT HLA scoring tsv
   - `phlat_qnsorted_bam`: Query name sorted bam of HLA hits only
   - `phlat_hla_map`: Newline separated list of HLA contigs mapped to

requirements:
- class: SubworkflowFeatureRequirement
inputs:
  # common
  output_basename: {type: string, doc: "String to use as basename for output file"}
  orientation: {type: ['null', {type: enum, name: orientation, symbols: ["--fr", "--rf", "--ff"]}], default: "--fr", doc: "If paired-end,
      read orientation"}
  maxins: {type: 'int?', doc: "maximum fragment length", default: 450}
  no_mixed: {type: 'boolean?', doc: "suppress unpaired alignments for paired reads", default: true}
  no_discordant: {type: 'boolean?', doc: "suppress discordant alignments for paired reads", default: true}
  no_unal: {type: 'boolean?', doc: "suppress unaligned reads", default: true}
  bwt2_threads: {type: 'int?', doc: "Num processing threads to use for bowtie2", default: 8}
  bwt2_ram: {type: 'int?', doc: "Num GB memory to make available for bowtie2", default: 16}
  extend_fail: {type: 'int?', doc: "give up extending after <int> failed extends in a row", default: 20}
  repetitive_seeds: {type: 'int?', doc: "for reads w/ repetitive seeds, try <int> sets of seeds", default: 3}
  max_seed_mm: {type: 'int?', doc: "max # mismatches in seed alignment; can be 0 or 1", default: 0}
  seed_substr_len: {type: 'int?', doc: "length of seed substrings; must be >3, <32", default: 20}
  seed_intervals: {type: 'string?', doc: "interval between seed substrings w/r/t read len (S,1,1.15)", default: "S,1,0.50"}
  bwt2_index: {type: File, doc: "bowtie2 tar-gzipped index"}
  bwt2_index_path_prefix: {type: string, doc: "unzipped path + prefix of index, i.e. index4phlat/ucsc.artHLA"}
  unpaired_reads: {type: 'File?', doc: "fastq reads if NOT paired"}
  reads: {type: 'File?', doc: "read1 fastq reads if paired-end"}
  mates: {type: 'File?', doc: "mates (read2) fastq reads if paired-end"}
  bwt2_header: {type: File, doc: "Desired header file to prepend to processed sorted bam"}
  preload_pickle: {type: File, doc: "pickle file from PHLAT resources"}
outputs:
  phlat_hla_summary: {type: File, doc: "Final PHLAT HLA scoring tsv", outputSource: phlat_hla_scoring/phlat_HLA_summary}
  phlat_qnsorted_bam: {type: File, doc: "Query name sorted bam of HLA hits only", outputSource: phlat_format_sam/phlat_name_sort}
  phlat_hla_map: {type: File, doc: "Newline separated list of HLA contigs mapped to", outputSource: phlat_format_sam/hla_map}
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
sbg:license: Apache License 2.0
sbg:publisher: KFDRC
