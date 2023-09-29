#$ -S /usr/bin/sh
#$ -j y
#$ -V
#$ -l m_mem_free=32G,h_vmem=80G
#$ -pe smp 2
#$ -cwd

module load python/2.7
module load sam-bcf-tools

phlatdir=/mnt/isilon/maris_lab/target_nbl_ngs/PHLAT/phlat-1.0_orig
datadir=/mnt/isilon/maris_lab/target_nbl_ngs/PHLAT/phlat-1.0_orig/example
indexdir=/mnt/isilon/maris_lab/target_nbl_ngs/PHLAT/phlat-1.0_orig/index4phlat
rsdir=/mnt/isilon/maris_lab/target_nbl_ngs/PHLAT/phlat-1.0_orig/results
b2url=/mnt/isilon/maris_lab/target_nbl_ngs/FormerLabMembers/Olivia/install/bowtie2-2.2.6/bowtie2
tag="example"
fastq1=${tag}"_1.fastq.gz"
fastq2=${tag}"_2.fastq.gz"
num_threads=2


python2.7 ${phlatdir}/dist/PHLAT.py -1 ${datadir}/${fastq1} -2 ${datadir}/${fastq2} -index $indexdir -b2url $b2url -tag $tag -e $phlatdir -o $rsdir -p $num_threads


