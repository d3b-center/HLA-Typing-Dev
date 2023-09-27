phlatdir=~/phlat-1.0
datadir=~/phlat-1.0/example
indexdir=~/phlat-1.0/index4phlat
rsdir=~/phlat-1.0/results
b2url=/Data/Programs/bin/bowtie2
tag="example"
fastq1=${tag}"_1.fastq.gz"
fastq2=${tag}"_2.fastq.gz"

python2.7 -O ${phlatdir}/dist/PHLAT.py -1 ${datadir}/${fastq1} -2 ${datadir}/${fastq2} -index $indexdir -b2url $b2url -tag $tag -e $phlatdir -o $rsdir


