#!/usr/bin/python2
import pysam
from extensions import *
import cPickle as pickle, os

def getAllSAMrecords(bamfile):
    reader = pysam.Samfile(bamfile, 'rb')
    iterAll = reader.fetch(until_eof=True)
    SAMrecordList = []
    for read in iterAll:
        if len(SAMrecordList) > 0:
            if read.qname == SAMrecordList[(-1)][0].qname:
                SAMrecordList[(-1)].append(SAMrecord(read))
            else:
                SAMrecordList.append([SAMrecord(read)])
        else:
            SAMrecordList.append([SAMrecord(read)])

    reader.close()
    return SAMrecordList


def getOneSAMrecord(bamfile, whichone):
    reader = pysam.Samfile(bamfile, 'rb')
    out = []
    counter = 1
    for read in reader.fetch(until_eof=True):
        counter = counter + 1
        if isinstance(whichone, int) and counter == whichone:
            out.append(SAMrecord(read))
        elif isinstance(whichone, str) and read.qname == whichone:
            out.append(SAMrecord(read))

    reader.close()
    return out


def initReadMaster(bamfile, HLAdb, SNPdb):
    reader = pysam.Samfile(bamfile, 'rb')
    iterAll = reader.fetch(until_eof=True)
    readmaster = readMaster()
    OneSAMrecordList = []
    for read in iterAll:
        if len(OneSAMrecordList) > 0:
            if read.qname == OneSAMrecordList[0].qname:
                OneSAMrecordList.append(SAMrecord(read))
            else:
                readmaster.addReads(OneSAMrecordList, HLAdb, SNPdb, 0.99)
                OneSAMrecordList = [SAMrecord(read)]
            readmaster.addReads(OneSAMrecordList, HLAdb, SNPdb, 0.99)
            OneSAMrecordList = []
        else:
            OneSAMrecordList.append(SAMrecord(read))

    readmaster.setHeader(reader.header)
    reader.close()
    readmaster.removeExistingBadReads()
    return readmaster


def go(tag, outdir, phlatdir, rmdir=True, change=0):
    bamfile = outdir + '/' + tag + '.tmp/' + tag + '_HLA.qnsorted.bam'
    outfile = outdir + '/' + tag + '_HLA.sum'
    mapfile = outdir + '/' + tag + '.tmp/' + tag + '_HLA.map'
    f1 = file(phlatdir + '/dist/resources/preload1.pk', 'rb')
    hladb, snpdb = pickle.load(f1)
    f1.close()
    hladb.updateHLAdb(findAllelesToSearch(mapfile))
    snpdb.updateSNPdb(hladb)
    readmaster = initReadMaster(bamfile, hladb, snpdb)
    allelemaster = alleleMaster(readmaster.getReadDict(), hladb, mapfile)
    allelemaster.analysisByLevels(hladb, keepslim=True)
    hladb.updateHLAdb(allelemaster.getAllelesToSearch())
    snpdb.updateSNPdb(hladb)
    readmaster.retainReadsByName(allelemaster.getReadNamesToKeep())
    pileupper = readPileupper()
    pileupper.pileupBySample(readmaster.getReadDict(), snpdb)
    LLgenoDict = calLLgenoDict(pileupper.extbasemap, doRound=True)
    llmaster = LLmaster(LLgenoDict, pileupper.phasemap, allelemaster, hladb, snpdb)
    outputResults(outfile, llmaster)
    os.system('rm ' + outdir + '/' + tag + '.tmp' + '/*.map')
    if rmdir:
        os.system('rm -r ' + outdir + '/' + tag + '.tmp')