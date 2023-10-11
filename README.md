# HLA Typing Development

This repo is for testing and benchmarking various HLA typing software.
In it's current state, it is in alpha phase and should only be used if you are familiar with this and know what you are doing

## [PHLAT](docs/D3b_PHLAT_README.md)

This is a software revival of sorts based on the the work from [PHLAT: Inference of High-Resolution HLA Types from RNA and Whole Exome Sequencing](https://pubmed.ncbi.nlm.nih.gov/29858810/). While the link the authors provided no longer works, our collaborators have found this tool invaluable and have provided us their copy.
The `phlat-1.0` directory is there for a point of reference, please refer to the [D3b PHLAT README](docs/logo/D3b_PHLAT_README.md) for our functionally equivalent implementation


### Usage Params
```
Usage: PHLAT.py -1 fastq1 [-2 fastq2] -index indexdir -b2url b2url [-tag samplename] [-p nthreads] [-e phlatdir] [-o outdir] [-pe 1]
-1: fastq file of the reads if single-end, or the first reads if paired-end
-2: fastq file of the second reads if paired-end;ignored if single-end
-index: url to the index files for Bowtie2 [default: index4phlat subfolder in phlat-1.0 packge]
-b2url: url to Bowtie2 executable
-tag: name label for the sample associated with the fastq files
-p: number of threads for running Bowtie2 [default 8]
-e: url to the home folder of phlat-1.0 package
-o: url to a directory where results shall be stored
-pe: flag indicating whether the data shall be treated as paired-end(1) or single-end(0) [default 1]
```

## [T1K](docs/T1K_README.md)

