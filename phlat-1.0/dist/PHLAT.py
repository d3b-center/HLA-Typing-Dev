#!/usr/bin/python
import prepare
import utilities
from prepare import *
from utilities import *

import sys
import datetime
import os

############## Command line execution starts here ##################
if __name__=="__main__":
    if (len(sys.argv)<7):
       if len(sys.argv)!=2 or (len(sys.argv)==2 and sys.argv[1]!="-h"):
          print "ERROR: insufficient input parameters"
       print "Usage: PHLAT.py -1 fastq1 [-2 fastq2] -index indexdir -b2url b2url [-tag samplename] [-p nthreads] [-e phlatdir] [-o outdir] [-pe 1]"
       print "-1: fastq file of the reads if single-end, or the first reads if paired-end"
       print "-2: fastq file of the second reads if paired-end;ignored if single-end"
       print "-index: url to the index files for Bowtie2 [default: index4phlat subfolder in phlat-1.0 packge]"
       print "-b2url: url to Bowtie2 executable"
       print "-tag: name label for the sample associated with the fastq files"
       print "-p: number of threads for running Bowtie2 [default 8]"
       print "-e: url to the home folder of phlat-1.0 package"
       print "-o: url to a directory where results shall be stored"
       print "-pe: flag indicating whether the data shall be treated as paired-end(1) or single-end(0) [default 1]"
       sys.exit()
 
    args=sys.argv[1:]
    outdir=""
    phlatdir="../"
    tag="sample"
    nthreads=8
    artindex="index4phlat"
    b2url=""
    fastq1=""
    fastq2=""
    ispe=None

    while len(args)>0:
        switch=args.pop(0)
        val=args.pop(0)
        if switch=="-tag":
            tag=val
        elif switch=="-o":
            outdir=val
        elif switch=="-e":
            phlatdir=val
        elif switch=="-index":
            artindex=val
        elif switch=="-b2url":
            b2url=val
        elif switch=="-1":
            fastq1=val
        elif switch=="-2":
            fastq2=val
        elif switch=="-pe":
            ispe=int(val)
        elif switch=="-p":
            nthreads=int(val)
    prepare(fastq1,fastq2,artindex,b2url,tag,outdir,phlatdir,ispe,"--fr", nthreads)

    start=datetime.datetime.now()
    print "..... Running PHLAT ......."
    go(tag,outdir,phlatdir)
    end=datetime.datetime.now()
    print ".....Done! Total PHLAT process time:"+str(end-start)+".....\n"





