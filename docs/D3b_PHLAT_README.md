# D3b PHLAT Implementation

This is a software revival of sorts based on the the work from [PHLAT: Inference of High-Resolution HLA Types from RNA and Whole Exome Sequencing](https://pubmed.ncbi.nlm.nih.gov/29858810/).
It was unfortunately implemented in python 2.7 and currently there are no plans to attempt to port over the original code to python3.
However, the workflow implementation should make this far easier and more accessible to use.
Some of the support scripts, like `scripts/extensions.py` and `scripts/utilities.py` have been copied over as-is to help ensure functional equivalency, while the former `PHAT.py` (`go()` function replaced by `scripts/phlat_d3b.py`) and `prepare.py` have been replaced with cwltool implementations for ease of use and readability

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

