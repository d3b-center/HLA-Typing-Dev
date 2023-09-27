# uncompyle6 version 3.7.4
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.6.8 (default, Apr  2 2020, 13:34:55) 
# [GCC 4.8.5 20150623 (Red Hat 4.8.5-39)]
# Embedded file name: /root/phlat2/dist/prepare.py
# Compiled at: 2017-09-23 15:54:26
import sys, subprocess, os, pysam

def prepare(fastq1, fastq2, artindex, b2url, tag, outdir, phlatdir, ispe, orientation, nthreads=8):
    tmpdir = outdir + '/' + tag + '.tmp'
    srcdir = phlatdir + '/dist/resources'
    if fastq1 == '' and fastq2 == '':
        print 'ERROR: must provide at least one fastq file'
        print 'Usage: PHLAT.py -1 fastq1 [-2 fastq2] -index indexdir -b2url b2url [-orientation orientation] [-e phlatdir] [-o outdir] [-tag samplename] [-pe 1] [-rmRare 0]'
        sys.exit()
    if fastq1 == '' and fastq2 != '':
        fastq1 = fastq2
        fastq2 = ''
    if ispe is None:
        if fastq2 != '':
            ispe = 1
        else:
            ispe = 0
    if ispe and fastq2 == '':
        ispe = 0
        print 'WARNING: less than two fastq files provided for paired-end mode, will process as single-end'
    if not ispe and fastq1 != '' and fastq2 != '':
        print 'WARNING: two fastq files provided despite of single-end mode,will process both files by ignoring the paired constrains'
    try:
        if os.path.exists(tmpdir):
            subprocess.call('rm -r ' + tmpdir, shell=True)
        subprocess.call('mkdir ' + tmpdir, shell=True)
        print 'mkdir ' + tmpdir
        print '.....Process Bowtie 2 mapping on ' + tag + '.......\n'
        if ispe and fastq1 != '' and fastq2 != '':
            cmd = b2url + ' ' + str(orientation) + ' --maxins 450 --no-mixed --no-discordant --no-unal -p ' + str(nthreads) + ' -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x ' + artindex + '/ucsc.artHLA -1 ' + fastq1 + ' -2 ' + fastq2 + ' |grep -P "\t(A\\*|B\\*|C\\*|DQA1\\*|DQB1\\*|DRB1\\*)" |grep -P -v "\tchr" >' + tmpdir + '/' + tag + '.sam'
            subprocess.call(cmd, shell=True)
            print cmd
        if not ispe and fastq2 == '':
            cmd = b2url + ' ' + str(orientation) + ' --maxins 450 --no-mixed --no-discordant --no-unal -p ' + str(nthreads) + ' -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x ' + artindex + '/ucsc.artHLA -U ' + fastq1 + '|grep -P "\t(A\\*|B\\*|C\\*|DQA1\\*|DQB1\\*|DRB1\\*)" |grep -P -v "\tchr" > ' + tmpdir + '/' + tag + '.sam'
            subprocess.call(cmd, shell=True)
            print cmd
        if not ispe and fastq1 != '' and fastq2 != '':
            cmd = b2url + ' ' + str(orientation) + ' --maxins 450 --no-mixed --no-discordant --no-unal -p ' + str(nthreads) + ' -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x ' + artindex + '/ucsc.artHLA -U ' + fastq1 + '|grep -P "\t(A\\*|B\\*|C\\*|DQA1\\*|DQB1\\*|DRB1\\*)" |grep -P -v "\tchr" |awk  -F"\t" \'BEGIN {OFS="\t"} {$1=$1."_1"}{print $0}\' - >' + tmpdir + '/' + tag + '.tmp1.sam'
            subprocess.call(cmd, shell=True)
            cmd = b2url + ' ' + str(orientation) + ' --maxins 450 --no-mixed --no-discordant --no-unal -p ' + str(nthreads) + ' -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x ' + artindex + '/ucsc.artHLA -U ' + fastq2 + '|grep -P "\t(A\\*|B\\*|C\\*|DQA1\\*|DQB1\\*|DRB1\\*)" |grep -P -v "\tchr" |awk  -F"\t" \'BEGIN {OFS="\t"} {$1=$1."_2"}{print $0}\' - > ' + tmpdir + '/' + tag + '.tmp2.sam'
            subprocess.call(cmd, shell=True)
            print cmd
            subprocess.call('cat ' + tmpdir + '/' + tag + '.tmp1.sam ' + tmpdir + '/' + tag + '.tmp2.sam >' + tmpdir + '/' + tag + '.sam ', shell=True)
        print '.....Prepare files of ' + tag + ' for PHLAT.......\n'
        cmd = 'cut -f 3 ' + tmpdir + '/' + tag + '.sam |grep -P "(A|B|C|DQA1|DQB1|DRB1)" > ' + tmpdir + '/' + tag + '_HLA.map'
        subprocess.call(cmd, shell=True)
        print cmd
        cmd = 'awk  -F"\t" \'BEGIN {OFS="\t"} {if($3~/^A/){$4=$4-1+29910331;$8=$8-1+29910331}} {if($3~/^B/){$4=$4-1+31322260;$8=$8-1+31322260}}{if($3~/^C/){$4=$4-1+31236946;$8=$8-1+31236946}} {if($3~/^DQA1/){$4=$4-1+32605236;$8=$8-1+32605236}} {if($3~/^DQB1/){$4=$4-1+32628013;$8=$8-1+32628013}} {if($3~/^DRB1/){$4=$4-1+32546868;$8=$8-1+32546868}} {if($3~/^(A|B|C|D)/){$3="chr6"}} {$5="255"}{print $0}\' ' + tmpdir + '/' + tag + '.sam' + '|cat ' + srcdir + '/bowtie.header - >' + tmpdir + '/' + tag + '.tmp.sam'
        subprocess.call(cmd, shell=True)
        print cmd
        pysam.view('-S', tmpdir + '/' + tag + '.tmp.sam', '-b', '-o' + tmpdir + '/' + tag + '.tmp.bam')
        pysam.sort('-n', tmpdir + '/' + tag + '.tmp.bam', tmpdir + '/' + tag + '_HLA.qnsorted')
        #subprocess.call('rm ' + tmpdir + '/' + tag + '.tmp*', shell=True)
        #subprocess.call('rm ' + tmpdir + '/' + tag + '.sam*', shell=True)
    except:
        print 'execution failed\n'
        subprocess.call('rm -r ' + tmpdir, shell=True)

    return
