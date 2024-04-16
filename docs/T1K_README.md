# T1K Genotyper
[The ONE genotyper for Kir and HLA](https://github.com/mourisl/T1K/).
Currently implemented as a tool, the [usage guide](https://github.com/mourisl/T1K/tree/v1.0.5#usage) explains all the possible tool options.

## Inputs
### Reference
 - `t1k_ref_tar`: Tar ball with T1K reference
 - `t1k_ref_path`: String of path to un-tarred hla ref, i.e. kiridx/kiridx_rna_seq.fa
 - `output_basename`: Prefix string for output file names. Default inferred from input

 ### Input Reads

The reads for T1K can be delivered as either BAM or fastq. The [T1K documentation](https://github.com/mourisl/T1K/tree/v1.0.5#inputoutput) outlines the various ways in which these files can be delivered. These match this tool's inputs: 

 - `reads`: Port for `-1` input, read1 fastq reads if paired-end (omit if single end or interleaved)
 - `mates`: Port for `-2` input, mates (read2) fastq reads if paired-end (omit if single end or interleaved)
 - `single_end`: Port for `-u` input, fastq reads if single-end (omit if paired-end or interleaved)
 - `interleaved`: Port for `-i` input, fastq reads if interleaved (omit if single end, or read1 and mates are separate files)
 - `bam`: Port for `-b` input, if bam input (use if input reads are bam and not fastq)

### Presets

T1K can be used for RNA-seq, WES, or WGS analysis. The default parameters are for RNA-seq. If you wish to perform a different analysis, you can change the `preset` field. That field has the following options:
- `hla`: HLA genotyping in general 
- `wgs`: HLA genotyping on WGS data
- `kir-wgs`: KIR genotyping on WGS data
- `kir-wes`: KIR genotyping on WES data

For more information, see the [T1K documentation](https://github.com/mourisl/T1K/tree/v1.0.5#preset-parameters).

### Tooling Options

All tool options are available in this tool. The options are too numerous to detail here. See the T1K documentation for information when and how to use these parameters.

### Resource Management
 - `threads`: Num processing threads to use
 - `ram`: GB of memory to make available

## Outputs
 - `aligned_fasta`
 - `allele_tsv`
 - `allele_vcf`
 - `candidate_fastqs`
 - `genotype_tsv`
 - `read_assignments`

Please see the [T1K documentation](https://github.com/mourisl/T1K/tree/v1.0.5#inputoutput) for descriptions of the outputs.