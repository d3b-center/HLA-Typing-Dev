# T1K Genotyper
["The ONE genotyper for Kir and HLA](https://github.com/mourisl/T1K/).
Currently implemented as a tool, the [usage guide](https://github.com/mourisl/T1K/tree/v1.0.2#usage) explains all the possible tool option.s Currently, the cwltool is limited, but will be expanded with future development

## tools/t1k.cwl
### Required Inputs
 Input reads:
 - `reads`: read1 fastq reads if paired-end (omit if single end or interleaved)
 - `mates`: mates (read2) fastq reads if paired-end (omit if single end or interleaved)
 - `single_end`: fastq reads if single-end (omit if paired-end or interleaved)
 - `interleaved`: fastq reads if interleaved (omit if single end, or read1 and mates are separate files)
 - `bam`: if bam input (use if input reads are bam and not fastq)
 Reference:
 - `t1k_ref_tar`: Tar ball with T1K reference
 - `t1k_ref_path`: String of path to un-tarred hla ref, i.e. kiridx/kiridx_rna_seq.fa
 - `output_basename`: Prefix string for output file names. Default inferred from input
### Optional resource
 - `threads`: Num processing threads to use
 - `ram`: Num GB memory to make available

### Outputs
Described [here](https://github.com/mourisl/T1K/tree/v1.0.2#inputoutput)
 - `aligned_fasta`
 - `allele_tsv`
 - `allele_vcf`
 - `candidate_fastqs`
 - `genotype_tsv`