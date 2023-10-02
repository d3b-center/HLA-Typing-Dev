"""
A rewrite of the PHLAT.py go function to perform final scoring
"""
import argparse
from utilities import *
from extensions import *
import cPickle as pickle
import pdb

parser = argparse.ArgumentParser(
    description="Run HLA typing using PHLAT algorithm on preprocessed inputs"
)
parser.add_argument(
    "-q",
    "--qsorted-bam",
    action="store",
    dest="qsort_bam",
    help="Queryname sorted bam",
)
parser.add_argument(
    "-m", "--map-file", action="store", dest="hla_map", help="HLA map file"
)
parser.add_argument(
    "-p", "--preload-pickle", action="store", dest="preload_pickle", help="pickle file from PHLAT resources"
)
parser.add_argument(
    "-o",
    "--output-basename",
    action="store",
    dest="output_basename",
    help="output_basename",
)


args = parser.parse_args()

bamfile = args.qsort_bam
outfile = args.output_basename + '_HLA.sum'
mapfile = args.hla_map
f1 = file(args.preload_pickle, 'rb')
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
pdb.set_trace()
llmaster = LLmaster(LLgenoDict, pileupper.phasemap, allelemaster, hladb, snpdb)
outputResults(outfile, llmaster)
