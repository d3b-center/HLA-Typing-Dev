#!/usr/bin/python2
# uncompyle6 version 3.7.4
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.6.8 (default, Apr  2 2020, 13:34:55) 
# [GCC 4.8.5 20150623 (Red Hat 4.8.5-39)]
# Embedded file name: /home/yu.bai/phlat-release/dist/extensions.py
# Compiled at: 2014-07-07 12:23:59
import sets, copy, math

def convertName(aname, resolution='4d', verbal=True):
    resolution = resolution.lower()
    s = aname.split('*')
    a = s[0]
    if not verbal and '_' in s[0]:
        tmp = s[0].split('_')
        a = tmp[1]
    if len(s) <= 1:
        return a + '*'
    shortallelename = a + '*' + s[1]
    if ':' in s[1]:
        myd = s[1].split(':')
        if resolution == '2d':
            shortallelename = a + '*' + myd[0]
        if resolution == '4d':
            shortallelename = a + '*' + myd[0] + ':' + myd[1]
        if resolution == '6d' and len(myd) >= 3:
            shortallelename = a + '*' + myd[0] + ':' + myd[1] + ':' + myd[2]
        if resolution == '8d' and len(myd) >= 4:
            shortallelename = a + '*' + myd[0] + ':' + myd[1] + ':' + myd[2] + ':' + myd[3]
    if resolution == 'locus':
        shortallelename = a
    return shortallelename


def findall(s, ch, offset):
    return [ i + offset for i, ltr in enumerate(s) if ltr == ch ]


def gauss2pval(val, mean, var, counts):
    if len(counts) == 1 and counts[0] > 20:
        return 0.1
    else:
        if mean is None or var is None:
            return 0
        if len(counts) < 3:
            return 0.1
        zscore = 1.0 * (val - mean) / math.sqrt(var)
        return math.erfc(zscore / math.sqrt(2)) / 2


def findWhich(list, query, op):
    length = len(list)
    rs = []
    if op == '<=':
        for r in range(length):
            if list[r] <= query:
                rs.append(r)

    elif op == '<':
        for r in range(length):
            if list[r] < query:
                rs.append(r)

    elif op == '>=':
        for r in range(length):
            if list[r] >= query:
                rs.append(r)

    elif op == '>':
        for r in range(length):
            if list[r] > query:
                rs.append(r)

    return rs


def mymean(num):
    return sum(num) * 1.0 / len(num)


def myvar(num):
    if len(num) < 2:
        if len(num) < 1:
            return None
        else:
            return num[0]

    mu = mymean(num)
    return sum([ math.pow(x - mu, 2) for x in num ]) * 1.0 / (len(num) - 1)


def findNorm(counts):
    if len(counts) < 3 or max(counts) < 10:
        if len(counts) < 1 or max(counts) < 10:
            return [None, None]
        if len(counts) < 3:
            return [1, 1]
    counts = sorted(counts, reverse=True)
    maxval = counts[0]
    minval = counts[(-1)]
    if maxval == minval:
        return [1, 1]
    else:
        bins = [ x / 25.0 for x in range(minval * 25, maxval * 25 + 1, maxval - minval) ]
        hists = [
         0] * len(bins)
        for val in counts:
            y = max(findWhich(bins, val, '<='))
            hists[y] = hists[y] + 1

        smoothList(hists)
        smoothList(hists)
        peaks = []
        for i in range(len(hists)):
            ispeak = False
            if i > 0 and i < len(hists) - 1:
                if hists[i] >= hists[(i + 1)] and hists[i] >= hists[(i - 1)] and i >= (len(bins) - 1) * 0.5:
                    ispeak = True
            if ispeak:
                peaks.append(bins[i])

        if len(peaks) == 0:
            peaks = maxval
            meanMax = peaks
        else:
            meanMax = max(peaks)
        peakMax = bins.index(meanMax)
        lowhist = max(min(hists[len(hists) / 2:]), hists[peakMax] / 7.39)
        span = max(maxval - meanMax, meanMax - bins[max(findWhich(hists, lowhist, '<='))])
        rangeMax = [meanMax - span, maxval]
        out = []
        for i in counts:
            if i <= rangeMax[1] and i >= rangeMax[0]:
                out.append(i)

        return [
         mymean(out), max(1, myvar(out))]


def smoothList(inlist):
    for index in range(len(inlist)):
        if index > 0 and index < len(inlist) - 1:
            inlist[index] = (inlist[(index - 1)] + inlist[index] + inlist[(index + 1)]) / 3.0
        elif index == 0:
            inlist[index] = (inlist[index] + inlist[(index + 1)]) / 2.0
        else:
            inlist[index] = (inlist[index] + inlist[(index - 1)]) / 2.0


def isValidBase(char):
    if char.lower() == 'a' or char.lower() == 't' or char.lower() == 'c' or char.lower() == 'g':
        return True
    return False


def readAllelesToSearchFile(AllelesToSearchFile):
    AllelesToSearch = sets.Set([])
    if AllelesToSearchFile is not None:
        filehandle = open(AllelesToSearchFile, 'r')
        for strline in filehandle:
            strline = strline.strip()
            safe = strline.split('\t')
            AllelesToSearch.add(safe[0])

        filehandle.close()
    return AllelesToSearch


def getClippedReadAlignmentPos(read):
    mystart = read.pos + 1
    tuplelist = read.cigar
    matchlen = 0
    for each in tuplelist:
        if each[0] == 0 or each[0] == 2 or each[0] == 3:
            matchlen = matchlen + each[1]

    myend = mystart + matchlen - 1
    return [mystart, myend]


def cigarReader(read):
    tuplelist = read.cigar
    formattedRead = ''
    seq = str(read.seq)
    currentChar = 0
    matchlen = 0
    for each in tuplelist:
        if each[0] == 0:
            formattedRead = formattedRead + seq[currentChar:currentChar + each[1]]
            currentChar = currentChar + each[1]
            matchlen = matchlen + each[1]
        elif each[0] == 1 or each[0] == 4:
            currentChar = currentChar + each[1]
        elif each[0] == 2 or each[0] == 3:
            matchlen = matchlen + each[1]
            deletion = ''
            for j in xrange(each[1]):
                deletion = deletion + 'D'

            formattedRead = formattedRead + deletion

    mystart = read.pos + 1
    myend = mystart + matchlen - 1
    return [formattedRead, mystart, myend]


class SAMrecord():

    def __init__(self, read):
        formattedRead = ''
        seq = str(read.seq)
        currentChar = 0
        formattedQual = []
        self.alignStart = long(read.pos + 1)
        self.alignEnd = read.aend
        self.DposInRead = sets.Set([])
        self.ntotal = 0
        self.nmatches = 0
        self.qname = read.qname
        self.seq = read.seq
        for each in read.cigar:
            if each[0] == 0:
                formattedRead = formattedRead + seq[currentChar:currentChar + each[1]]
                formattedQual.extend(map(lambda x: ord(x) - 33, read.qual[currentChar:currentChar + each[1]]))
                currentChar = currentChar + each[1]
            elif each[0] == 1 or each[0] == 4:
                currentChar = currentChar + each[1]
            elif each[0] == 2 or each[0] == 3:
                tmpsite = len(formattedRead)
                deletion = ''
                for j in xrange(each[1]):
                    deletion = deletion + 'D'

                formattedRead = formattedRead + deletion
                formattedQual.extend(map(lambda x: ord(x) - 33, read.qual[currentChar:currentChar + each[1]]))
                self.DposInRead = self.DposInRead.union(xrange(self.alignStart + tmpsite, self.alignStart + tmpsite + each[1]))

        self.formattedRead = formattedRead
        self.formattedQual = formattedQual
        self.whichMaxIdentityAlleles = []
        self.maxIdentity = 0.0

    def setFormattedRead(self, formattedRead):
        self.formattedRead = formattedRead

    def setAlignStart(self, alignStart):
        self.alignStart = alignStart

    def setAlignEnd(self, alignEnd):
        self.alignEnd = alignEnd

    def getFormattedRead(self):
        return self.formattedRead

    def getCharAt(self, pos):
        if pos > self.alignEnd or pos < self.alignStart:
            return
        newpos = pos - self.alignStart
        return self.formattedRead[newpos].upper()
        return

    def getQualAt(self, pos):
        if pos > self.alignEnd or pos < self.alignStart:
            return
        newpos = pos - self.alignStart
        return self.formattedQual[newpos]
        return

    def getAlignStart(self):
        return self.alignStart

    def getAlignEnd(self):
        return self.alignEnd

    def clearAlleleRelated(self):
        self.whichMaxIdentityAlleles = []
        self.maxIdentity = 0.0

    def setWhichMaxIdentityAlleles(self, whichalleles):
        self.whichMaxIdentityAlleles = whichalleles

    def clearWhichMaxIdentityAlleles(self):
        self.whichMaxIdentityAlleles = []

    def addWhichMaxIdentityAlleles(self, newwhichalleles):
        self.whichMaxIdentityAlleles = list(sets.Set(self.whichMaxIdentityAlleles).union(sets.Set(newwhichalleles)))

    def getWhichMaxIdentityAlleles(self):
        return self.whichMaxIdentityAlleles

    def setMaxIdentity(self, maxidentity):
        self.maxIdentity = float(maxidentity)

    def getMaxIdentity(self):
        return self.maxIdentity


def calMatches(thisread, thisallele, readstart1, readstart2, allelestart, mattersnp1, mattersnp2):
    nmatches = 0
    ntotal = 0
    for pos in mattersnp1:
        ntotal = ntotal + 1
        if thisread[0][(pos - readstart1)] == thisallele[(pos - allelestart)]:
            nmatches = nmatches + 1
        if ntotal - nmatches > MAXMIS:
            break

    if len(thisread) > 1 and ntotal - nmatches <= MAXMIS:
        for pos in mattersnp2:
            ntotal = ntotal + 1
            if thisread[1][(pos - readstart2)] == thisallele[(pos - allelestart)]:
                nmatches = nmatches + 1
            if ntotal - nmatches > MAXMIS:
                break

    return float(nmatches) / ntotal


class HLAdatabase():

    def __init__(self, HLAseqFile, AlleleFreqFile, mapfile=None):
        self.HLAseqs = []
        self.starts = []
        self.stops = []
        self.names = []
        self.minstartpos = float('inf')
        self.maxstoppos = -1
        self.Locus2Indices = {}
        self.Dpos = sets.Set([])
        self.MaxAlleleFreqs = {}
        self.setAlleleFreqs(AlleleFreqFile)
        AllelesToSearch = findAllelesToSearch(mapfile)
        self.setHLAdbByAllelesToSearch(HLAseqFile, AllelesToSearch)
        self.recordDpositions()

    def clear(self):
        self.HLAseqs = []
        self.starts = []
        self.stops = []
        self.names = []
        self.minstartpos = float('inf')
        self.maxstoppos = -1
        self.Locus2Indices = {}
        self.Dpos = sets.Set([])
        self.MaxAlleleFreqs = {}

    def setHLAdbByAllelesToSearch(self, HLAseqFile, rmRare, AllelesToSearch):
        filehandle = open(HLAseqFile, 'r')
        counter = 0
        for strline in filehandle:
            strline = strline.strip()
            s = strline.split('\t')
            tmp = convertName(s[0], '', verbal=False)
            if len(AllelesToSearch) == 0 or len(AllelesToSearch) > 0 and (tmp in AllelesToSearch or s[0] in AllelesToSearch):
                locus = convertName(s[0], 'locus')
                if locus not in self.Locus2Indices:
                    self.Locus2Indices[locus] = []
                self.Locus2Indices[locus].append(counter)
                self.HLAseqs.append(s[3])
                self.names.append(s[0])
                self.starts.append(long(s[1]))
                self.stops.append(long(s[2]))
                counter = counter + 1

        self.minstartpos = min(self.starts)
        self.maxstoppos = max(self.stops)
        filehandle.close()

    def updateHLAdb(self, AllelesToSearch, rmRare=0):
        rmindices = []
        for which in xrange(len(self.names)):
            if rmRare and self.MaxAlleleFreqs[convertName(self.names[which])] <= 0.0001 or len(AllelesToSearch) > 0 and convertName(self.names[which], '', False) not in AllelesToSearch and convertName(self.names[which], '') not in AllelesToSearch:
                rmindices.append(which)

        rmindices = sorted(rmindices, reverse=True)
        for qq in rmindices:
            del self.HLAseqs[qq]
            del self.starts[qq]
            del self.stops[qq]
            del self.names[qq]

        self.minstartpos = min(self.starts)
        self.maxstoppos = max(self.stops)
        self.Locus2Indices = {}
        for which in xrange(len(self.names)):
            locus = convertName(self.names[which], 'locus')
            if locus not in self.Locus2Indices:
                self.Locus2Indices[locus] = []
            self.Locus2Indices[locus].append(which)

        self.recordDpositions()

    def getHLAseqs(self):
        return self.HLAseqs

    def getNames(self):
        return self.names

    def getStartPositions(self):
        return self.starts

    def getStopPositions(self):
        return self.stops

    def getLocus2Indices(self):
        return self.Locus2Indices

    def getIndex2SeqByLocus(self, readstart, readend):
        index2seq = {}
        for index in self.getIndicesByLocus(readstart, readend):
            index2seq[index] = self.HLAseqs[index]

        return index2seq

    def getIndicesByLocus(self, readstart, readend):
        if not (readend < 29910331 or readstart > 29913232) and 'HLA_A' in self.Locus2Indices:
            return self.Locus2Indices['HLA_A']
        if not (readend < 31322260 or readstart > 31324935) and 'HLA_B' in self.Locus2Indices:
            return self.Locus2Indices['HLA_B']
        if not (readend < 31236946 or readstart > 31239848) and 'HLA_C' in self.Locus2Indices:
            return self.Locus2Indices['HLA_C']
        if not (readend < 32605236 or readstart > 32610541) and 'HLA_DQA1' in self.Locus2Indices:
            return self.Locus2Indices['HLA_DQA1']
        if not (readend < 32628013 or readstart > 32634384) and 'HLA_DQB1' in self.Locus2Indices:
            return self.Locus2Indices['HLA_DQB1']
        if not (readend < 32407728 or readstart > 32411687) and 'HLA_DRA' in self.Locus2Indices:
            return self.Locus2Indices['HLA_DRA']
        else:
            if not (readend < 32546868 or readstart > 32557519) and 'HLA_DRB1' in self.Locus2Indices:
                return self.Locus2Indices['HLA_DRB1']
            if not (readend < 30457309 or readstart > 30460358) and 'HLA_E' in self.Locus2Indices:
                return self.Locus2Indices['HLA_E']
            if not (readend < 29691241 or readstart > 29694184) and 'HLA_F' in self.Locus2Indices:
                return self.Locus2Indices['HLA_F']
            if not (readend < 33036427 or readstart > 33041347) and 'HLA_DPA1' in self.Locus2Indices:
                return self.Locus2Indices['HLA_DPA1']
            if not (readend < 33043819 or readstart > 33054015) and 'HLA_DPB1' in self.Locus2Indices:
                return self.Locus2Indices['HLA_DPB1']
            if not (readend < 32916641 or readstart > 32920813) and 'HLA_DMA' in self.Locus2Indices:
                return self.Locus2Indices['HLA_DMA']
            if not (readend < 32902748 or readstart > 32908584) and 'HLA_DMB' in self.Locus2Indices:
                return self.Locus2Indices['HLA_DMB']
            if not (readend < 32974615 or readstart > 32977313) and 'HLA_DOA' in self.Locus2Indices:
                return self.Locus2Indices['HLA_DOA']
            if not (readend < 32780993 or readstart > 32784728) and 'HLA_DOB' in self.Locus2Indices:
                return self.Locus2Indices['HLA_DOB']
            return []

    def getIndexByAlleleName(self, aname):
        if aname in self.names:
            return self.names.index(aname)
        else:
            return -1

    def recordDpositions(self):
        Dpos = sets.Set([])
        for each in xrange(len(self.HLAseqs)):
            Dpos = sets.Set(Dpos).union(findall(self.HLAseqs[each], 'D', self.starts[each]))

        self.Dpos = Dpos

    def getDpositions(self):
        return self.Dpos

    def setAlleleFreqs(self, AlleleFreqFile):
        filehandle = open(AlleleFreqFile, 'r')
        counter = 0
        for strline in filehandle:
            counter = counter + 1
            strline = strline.strip()
            s = strline.split('\t')
            if counter == 1:
                for i in xrange(len(s) - 1):
                    pass

            else:
                for i in xrange(len(s) - 1):
                    if 'HLA_' + s[0] not in self.MaxAlleleFreqs:
                        self.MaxAlleleFreqs['HLA_' + s[0]] = float(s[(i + 1)])
                    elif float(self.MaxAlleleFreqs[('HLA_' + s[0])]) < float(s[(i + 1)]):
                        self.MaxAlleleFreqs['HLA_' + s[0]] = float(s[(i + 1)])

        filehandle.close()

    def getMaxFreqs(self):
        return self.MaxAlleleFreqs


class SNPdatabase():

    def __init__(self, SNPsiteFile, HLAdb):
        self.SNPsites = sets.Set([])
        self.SNPsitesNoD = sets.Set([])
        self.setSNPsites(SNPsiteFile, HLAdb)
        self.regions = initCoverRegion('CDS')
        self.withinRegionSites = sets.Set([])
        for qq in self.regions:
            self.withinRegionSites = self.withinRegionSites.union(sets.Set(xrange(qq[0], qq[1] + 1)))

    def updateSNPdb(self, HLAdb):
        self.SNPsitesNoD = sets.Set(self.SNPsites).difference(HLAdb.getDpositions())

    def setSNPsites(self, SNPsiteFile, HLAdb):
        Dpos = sets.Set([])
        if HLAdb is not None:
            Dpos = HLAdb.getDpositions()
        snpfilehandle = open(SNPsiteFile, 'r')
        for strline in snpfilehandle:
            strline = strline.strip()
            s = strline.split('\t')
            pos = long(s[0])
            indices = HLAdb.getIndicesByLocus(pos, pos)
            if len(indices) > 0 and float(s[8]) > 0:
                self.SNPsites.add(pos)
                if pos not in Dpos:
                    self.SNPsitesNoD.add(pos)

        snpfilehandle.close()
        return

    def getSNPsites(self):
        return self.SNPsites

    def getSNPsitesNoD(self):
        return self.SNPsitesNoD

    def setCoverRegion(self, mytype):
        self.regions = initCoverRegion(mytype)

    def getCoverRegion(self):
        return self.regions

    def isWithinRegion(self, pos):
        return pos in self.withinRegionSites

    def getSNPsitesWithin(self, minpos, maxpos):
        return sets.Set(xrange(minpos, maxpos + 1)).intersection(self.SNPsites)

    def getSNPsitesNoDWithin(self, minpos, maxpos):
        return sets.Set(xrange(minpos, maxpos + 1)).intersection(self.SNPsitesNoD)

    def getDnumByPos(self, pos, indices, starts, seqs):
        numD = 0
        start = starts[indices[0]]
        newpos = pos - start
        for each in indices:
            sq = seqs[each]
            if sq[newpos] == 'D':
                numD = numD + 1
                break

        return numD


class readMaster():

    def __init__(self):
        self.readDict = {}
        self.currentReadName = ''
        self.badreadnames = sets.Set([])
        self.header = {}

    def setHeader(self, header):
        self.header = header

    def getCurrentReadName(self):
        return self.currentReadName

    def getReadDict(self):
        return self.readDict

    def recordMaxIdentityAll(self, HLAdb, SNPdb):
        for mykey in self.readDict.keys():
            recordMaxIdentityOne(self.readDict[mykey], HLAdb, SNPdb)

    def setBadReadsByIdentity(self, minIdentity):
        self.badreadnames = sets.Set([])
        for mykey in self.readDict.keys():
            if self.readDict[mykey][0].getMaxIdentity() < minIdentity:
                self.badreadnames.add(mykey)

    def removeExistingBadReads(self):
        self.removeReadsByName(self.badreadnames)

    def removeReadsByName(self, readnames):
        for mykey in readnames:
            if mykey in self.readDict.keys():
                del self.readDict[mykey]

    def removeReads(self, reads):
        for read in reads:
            mykey = read.qname
            if mykey in self.readDict.keys():
                del self.readDict[mykey]

    def retainReadsByName(self, readnames):
        rmnames = sets.Set(self.readDict.keys()).difference(sets.Set(readnames))
        self.removeReadsByName(rmnames)

    def addReads(self, OneSAMrecordList, HLAdb, SNPdb, minIdentity, discardbad=True):
        if len(OneSAMrecordList) > 1:
            OneSAMrecordList = sorted(OneSAMrecordList, key=lambda check: check.alignStart)
        readname = OneSAMrecordList[0].qname
        recordMaxIdentityOne(OneSAMrecordList, HLAdb, SNPdb)
        if OneSAMrecordList[0].getMaxIdentity() < minIdentity:
            if not discardbad:
                self.currentReadName = readname
                self.badreadnames.add(readname)
                self.readDict[readname] = OneSAMrecordList
        else:
            self.currentReadName = readname
            self.readDict[readname] = OneSAMrecordList


def recordMaxIdentityOne(OneSAMrecordList, HLAdb, SNPdb):
    MAXMIS = 3
    thisread = [OneSAMrecordList[0].getFormattedRead()]
    readstart1 = OneSAMrecordList[0].getAlignStart()
    readstop1 = OneSAMrecordList[0].getAlignEnd()
    if len(OneSAMrecordList) > 1:
        readstart2 = OneSAMrecordList[1].getAlignStart()
        readstop2 = OneSAMrecordList[1].getAlignEnd()
        thisread.append(OneSAMrecordList[1].getFormattedRead())
    else:
        readstart2 = readstart1
        readstop2 = readstop1
    maxIdentity = 0.0
    whichmaxalleles = []
    myappend = whichmaxalleles.append
    maxntotal = 0
    maxnmatches = 0
    readMin = min([readstart1, readstart2, readstop1, readstop2])
    readMax = max([readstart1, readstart2, readstop1, readstop2])
    index2seq = HLAdb.getIndex2SeqByLocus(readMin, readMax)
    allelestart = HLAdb.getStartPositions()[index2seq.keys()[0]]
    allelestop = HLAdb.getStopPositions()[index2seq.keys()[0]]
    mattersnp1 = SNPdb.getSNPsitesNoDWithin(max(allelestart, readstart1), min(readstop1, allelestop)).intersection(SNPdb.withinRegionSites).difference(OneSAMrecordList[0].DposInRead)
    mattersnp2 = sets.Set([])
    if len(OneSAMrecordList) > 1:
        mattersnp2 = SNPdb.getSNPsitesNoDWithin(max(allelestart, readstart2), min(readstop2, allelestop)).intersection(SNPdb.withinRegionSites).difference(OneSAMrecordList[1].DposInRead)
    if len(mattersnp1) + len(mattersnp2) > 0:
        for index in index2seq:
            thisallele = index2seq[index]
            nmatches = 0
            ntotal = 0
            for pos in mattersnp1:
                ntotal = ntotal + 1
                if thisread[0][(pos - readstart1)] == thisallele[(pos - allelestart)]:
                    nmatches = nmatches + 1
                if ntotal - nmatches > MAXMIS:
                    break

            if len(thisread) > 1 and ntotal - nmatches <= MAXMIS:
                for pos in mattersnp2:
                    ntotal = ntotal + 1
                    if thisread[1][(pos - readstart2)] == thisallele[(pos - allelestart)]:
                        nmatches = nmatches + 1
                    if ntotal - nmatches > MAXMIS:
                        break

            nidentity = float(nmatches) / ntotal
            if nidentity > maxIdentity:
                maxIdentity = nidentity
                maxntotal = ntotal
                maxnmatches = nmatches
                whichmaxalleles = [index]
            elif nidentity == maxIdentity:
                whichmaxalleles.append(index)

    OneSAMrecordList[0].setWhichMaxIdentityAlleles(whichmaxalleles)
    OneSAMrecordList[0].setMaxIdentity(maxIdentity)
    OneSAMrecordList[0].ntotal = maxntotal
    OneSAMrecordList[0].nmatches = maxnmatches


class alleleMaster():

    def __init__(self, readDict, HLAdb, mapfile):
        self.alleleCountfd = {}
        self.alleleCount4d = {}
        self.locusCountfd = {}
        self.locusCount4d = {}
        self.allele2Reads = {}
        self.locus2Allele4d = {}
        self.allelesToSearch = sets.Set([])
        self.pvals = {}
        self.readNamesToKeep = sets.Set([])
        self.readlength = 0
        self.maptable = {}
        self.setMaster(readDict, HLAdb, mapfile)
        self.measure = 1000

    def setMasterByCopy(self, cp):
        self.alleleCountfd = cp.alleleCountfd
        self.alleleCount4d = cp.alleleCount4d
        self.locusCountfd = cp.locusCountfd
        self.locusCount4d = cp.locusCount4d
        self.allele2Reads = cp.allele2Reads
        self.locus2Allele4d = cp.locus2Allele4d
        self.allelesToSearch = cp.allelesToSearch
        self.readNamesToKeep = cp.readNamesToKeep
        self.pvals = cp.pvals
        self.maptable = cp.maptable
        self.readlength = cp.readlength
        self.measure = cp.measure

    def getReadNamesToKeep(self):
        return self.readNamesToKeep

    def getAllelesToSearch(self):
        return self.allelesToSearch

    def getAlleleCountfd(self):
        return self.alleleCountfd

    def getAlleleCount4d(self):
        return self.alleleCount4d

    def getLocusCount(self):
        return self.locusCount4d

    def clearAllele2Reads(self):
        self.allele2Reads = {}

    def setMaster(self, readDict, HLAdb, mapfile):
        self.maptable = readMapFile(mapfile)
        myget = readDict.get
        tmpread = readDict[readDict.keys()[0]][0]
        self.readlength = len(tmpread.seq)
        for onekey in readDict:
            oneread = readDict[onekey]
            readname = oneread[0].qname
            self.readNamesToKeep.add(readname)
            fourDigitAlleles = sets.Set([])
            whichalleles = oneread[0].getWhichMaxIdentityAlleles()
            for alleleindex in whichalleles:
                allele = HLAdb.getNames()[alleleindex]
                locus = convertName(allele, resolution='locus')
                allele4d = convertName(allele, resolution='4d')
                allele = convertName(allele, resolution='')
                self.alleleCountfd[allele] = self.alleleCountfd.get(allele, 0) + 1
                if locus not in self.locus2Allele4d:
                    self.locus2Allele4d[locus] = sets.Set([])
                if allele4d not in fourDigitAlleles:
                    fourDigitAlleles.add(allele4d)
                    self.locus2Allele4d[locus].add(allele4d)
                    if allele4d not in self.allele2Reads:
                        self.allele2Reads[allele4d] = sets.Set([])
                    self.allele2Reads[allele4d].add(readname)
                if self.alleleCountfd[allele] > 1 and allele not in self.allelesToSearch:
                    self.allelesToSearch.add(allele)

            if len(whichalleles) > 0:
                allele = HLAdb.getNames()[whichalleles[0]]
                locus = convertName(allele, resolution='locus')
                if locus not in self.locusCountfd:
                    self.locusCountfd[locus] = 1
                    self.locusCount4d[locus] = 1
                else:
                    self.locusCountfd[locus] = self.locusCountfd[locus] + 1
                    self.locusCount4d[locus] = self.locusCount4d[locus] + 1
            for uallele4d in fourDigitAlleles:
                self.alleleCount4d[uallele4d] = self.alleleCount4d.get(uallele4d, 0) + 1

    def analysisByLevels(self, hladb, keepslim=True):
        goodReadNames = sets.Set([])
        badReadNames = sets.Set([])
        newAllelesToSearch = []
        keptAllele4d = sets.Set([])
        AllMaxAlleles2 = []
        sumCounts = 0
        for eachL in self.locus2Allele4d:
            sumCounts = sumCounts + self.locusCount4d[eachL]
            MaxAlleles0 = []
            secondMaxAlleles0 = []
            thirdMaxAlleles0 = []
            MaxAlleles1 = []
            secondMaxAlleles1 = []
            MaxAlleles2 = []
            maxcount00 = 0
            secondmaxcount00 = 0
            secondmaxcount01 = 0
            thirdmaxcount00 = 0
            thirdmaxcount01 = 0
            maxcount11 = 0
            secondmaxcount11 = 0
            maxcount22 = 0
            tmpcount = {}
            level0 = {}
            level1 = {}
            level2 = {}
            tmpset = self.locus2Allele4d[eachL]
            for eacha in self.locus2Allele4d[eachL]:
                currentcount = len(self.allele2Reads[eacha])
                level0[eacha] = currentcount
                if currentcount > maxcount00:
                    MaxAlleles0 = [
                     eacha]
                    maxcount00 = currentcount
                elif currentcount == maxcount00:
                    MaxAlleles0.append(eacha)

            keptAllele4d = keptAllele4d.union(sets.Set(MaxAlleles0).difference(getRAllele0(MaxAlleles0, maxcount00, hladb, self.maptable)))
            for MaxAllele in MaxAlleles0:
                tmpcount[MaxAllele] = maxcount00

            tmpset = tmpset.difference(sets.Set(MaxAlleles0))
            for eacha in tmpset:
                currentcount = len(self.allele2Reads[eacha])
                if currentcount > secondmaxcount00:
                    secondMaxAlleles0 = [
                     eacha]
                    secondmaxcount00 = currentcount
                elif currentcount == secondmaxcount00:
                    secondMaxAlleles0.append(eacha)

            for eacha in tmpset:
                for MaxAllele in MaxAlleles0:
                    self.allele2Reads[eacha] = self.allele2Reads[eacha].difference(self.allele2Reads[MaxAllele])

                currentcount = len(self.allele2Reads[eacha])
                level1[eacha] = currentcount
                if currentcount > maxcount11:
                    MaxAlleles1 = [
                     eacha]
                    maxcount11 = currentcount
                elif currentcount == maxcount11:
                    MaxAlleles1.append(eacha)
                if eacha in secondMaxAlleles0:
                    if currentcount > secondmaxcount01:
                        secondmaxcount01 = currentcount
                self.alleleCount4d[eacha] = currentcount

            if len(MaxAlleles1) <= 0:
                pass
            else:
                for sa in MaxAlleles1:
                    for resta in tmpset:
                        if convertName(sa, '2d') == convertName(resta, '2d'):
                            currentcount = level0[resta]
                            if currentcount > thirdmaxcount00:
                                thirdMaxAlleles0 = [
                                 resta]
                                thirdmaxcount00 = currentcount
                            elif currentcount == thirdmaxcount00:
                                thirdMaxAlleles0.append(resta)

                for sa in tmpset:
                    currentcount = level1[sa]
                    if sa not in MaxAlleles1:
                        if currentcount > secondmaxcount11:
                            secondMaxAlleles1 = [
                             sa]
                            secondmaxcount11 = currentcount
                        elif currentcount == secondmaxcount11:
                            secondMaxAlleles1.append(sa)
                    if sa in thirdMaxAlleles0:
                        if currentcount > thirdmaxcount01:
                            thirdmaxcount01 = currentcount

                shared2d = False
                for secondMaxAllele in MaxAlleles1:
                    for MaxAllele in MaxAlleles0:
                        if convertName(secondMaxAllele, '2d') == convertName(MaxAllele, '2d'):
                            shared2d = True
                            break

                    if shared2d:
                        break

                secondshared2d = True
                for secondMaxAllele in MaxAlleles1:
                    for secondMaxAllele0 in secondMaxAlleles0:
                        if convertName(secondMaxAllele, '2d') != convertName(secondMaxAllele0, '2d'):
                            secondshared2d = False
                            break

                    if not secondshared2d:
                        break

                firstshared2d = True
                for MaxAllele in MaxAlleles0:
                    for secondMaxAllele0 in secondMaxAlleles0:
                        if convertName(MaxAllele, '2d') != convertName(secondMaxAllele0, '2d'):
                            firstshared2d = False
                            break

                    if not firstshared2d:
                        break

                secondshared2d1 = True
                for secondMaxAllele in MaxAlleles1:
                    for secondMaxAllele1 in secondMaxAlleles1:
                        if convertName(secondMaxAllele, '2d') != convertName(secondMaxAllele1, '2d'):
                            secondshared2d1 = False
                            break

                    if not secondshared2d1:
                        break

                ratio = 1.0 * maxcount00 / self.locusCount4d[eachL]
                secondratio = 1.0 * maxcount11 / self.locusCount4d[eachL]
                secondratio1 = 1.0 * secondmaxcount11 / self.locusCount4d[eachL]
                sratio01 = 1.0 * secondmaxcount01 / maxcount00
                keepcutoff = 0.01
                if sratio01 >= keepcutoff and ratio + secondratio > 0.75 or secondratio > 0.1 and secondshared2d and not firstshared2d:
                    keptAllele4d = keptAllele4d.union(sets.Set(secondMaxAlleles0).difference(getRAllele0(secondMaxAlleles0, secondmaxcount00, hladb, self.maptable)))
                hetcutoff = 0.05
                sratio11 = maxcount11 * 1.0 / maxcount00
                if ratio > 0.9 and not shared2d:
                    if sratio11 > hetcutoff:
                        keptAllele4d = keptAllele4d.union(sets.Set(MaxAlleles1).difference(getRAllele1(MaxAlleles1, maxcount11, hladb, self.maptable)))
                    for eacha in tmpset:
                        if eacha not in keptAllele4d:
                            self.alleleCount4d[eachL] = 0

                else:
                    for MaxAllele in MaxAlleles0:
                        for secondMaxAllele in sets.Set(MaxAlleles1).difference(getRAllele1(MaxAlleles1, maxcount11, hladb, self.maptable)):
                            max2d = convertName(MaxAllele, '2d')
                            secondmax2d = convertName(secondMaxAllele, '2d')
                            if max2d == secondmax2d and sratio11 >= hetcutoff or max2d != secondmax2d and ratio < 0.9 and sratio11 >= hetcutoff:
                                keptAllele4d.add(secondMaxAllele)
                            tmpcount[secondMaxAllele] = maxcount11

                    tmpset = tmpset.difference(sets.Set(MaxAlleles1))
                    for eacha in tmpset:
                        for secondMaxAllele in MaxAlleles1:
                            self.allele2Reads[eacha] = self.allele2Reads[eacha].difference(self.allele2Reads[secondMaxAllele])

                        currentcount = len(self.allele2Reads[eacha])
                        level2[eacha] = currentcount
                        if eacha not in keptAllele4d:
                            self.alleleCount4d[eacha] = currentcount
                        if currentcount > maxcount22:
                            MaxAlleles2 = [
                             eacha]
                            maxcount22 = currentcount
                        elif currentcount == maxcount22:
                            MaxAlleles2.append(eacha)
                        tmpcount[eacha] = currentcount

                    AllMaxAlleles2.extend(sets.Set(MaxAlleles2).difference(getRAllele1(MaxAlleles2, maxcount22, hladb, self.maptable)))
                saved = sets.Set(MaxAlleles0)
                saved = saved.union(sets.Set(secondMaxAlleles0))
                saved = saved.union(sets.Set(MaxAlleles1))
                saved = saved.union(sets.Set(MaxAlleles2))
                for aa in saved:
                    if aa not in self.pvals:
                        if aa in MaxAlleles0:
                            cp = copy.deepcopy(level0)
                            for bb in MaxAlleles0:
                                del cp[bb]

                            model = findNorm(list(cp.values()))
                            self.pvals[aa] = gauss2pval(level0[aa], model[0], model[1], list(cp.values()))
                        elif aa in MaxAlleles1:
                            cp = copy.deepcopy(level1)
                            for bb in MaxAlleles1:
                                del cp[bb]

                            model = findNorm(list(cp.values()))
                            self.pvals[aa] = gauss2pval(level1[aa], model[0], model[1], list(cp.values()))
                        elif aa in secondMaxAlleles0:
                            cp = copy.deepcopy(level0)
                            for bb in MaxAlleles0:
                                del cp[bb]

                            for bb in secondMaxAlleles0:
                                del cp[bb]

                            model = findNorm(list(cp.values()))
                            self.pvals[aa] = gauss2pval(level0[aa], model[0], model[1], list(cp.values()))
                        elif aa in MaxAlleles2:
                            cp = copy.deepcopy(level2)
                            for bb in MaxAlleles2:
                                del cp[bb]

                            model = findNorm(list(cp.values()))
                            self.pvals[aa] = gauss2pval(level2[aa], model[0], model[1], list(cp.values()))

        tmpcopy = copy.deepcopy(secondMaxAlleles0)
        for my2ndmaxallele in tmpcopy:
            if my2ndmaxallele in keptAllele4d and (self.alleleCount4d[my2ndmaxallele] < secondmaxcount01 or self.alleleCount4d[my2ndmaxallele] == 0):
                keptAllele4d.remove(my2ndmaxallele)

        self.measure = 1.0 * sumCounts * self.readlength / countBaseByRegions(initCoverRegion('CDS'))
        for myallele in self.allelesToSearch:
            myallele4d = convertName(myallele, '4d')
            passratio = self.alleleCount4d[myallele4d] * 1.0 / self.locusCount4d[convertName(myallele, 'locus')]
            passcount = self.alleleCount4d[myallele4d]
            if myallele4d in keptAllele4d or myallele4d in AllMaxAlleles2 and passratio > 0.1 and passcount >= 100:
                if myallele4d not in keptAllele4d:
                    pass
                newAllelesToSearch.append(myallele)
                goodReadNames = goodReadNames.union(self.allele2Reads[myallele4d])
            else:
                badReadNames = badReadNames.union(self.allele2Reads[myallele4d])

        if keepslim:
            self.allele2Reads = {}
        self.allelesToSearch = sets.Set(newAllelesToSearch)
        del newAllelesToSearch
        self.readNamesToKeep = goodReadNames


LOG10_3 = math.log10(3.0)
LOG10_PCR_ERR = math.log10(0.0001)
LOG10_PCR_ERR_PERSTATE = LOG10_PCR_ERR - LOG10_3
LOG10_1_MINUS_PCR_ERR = math.log10(1.0 - 0.0001)
PLOIDY = 2
P_err = 0.01
P_correct = 1 - P_err
L_err = math.log10(P_err / 15)
L_correct = math.log10(P_correct)

def log10PofObsbaseGivenRefbase(obsbase, obsqual, refbase):
    if obsbase.lower() == refbase.lower():
        logP = math.log10(1.0 - math.pow(10, obsqual / -10.0))
    else:
        logP = obsqual / -10.0 + -LOG10_3
    return logP


def log10FourRefbase(obsbase, obsqual):
    log10Pby4Refbase = {}
    log10Pby4Refbase['A'] = 0
    log10Pby4Refbase['C'] = 0
    log10Pby4Refbase['G'] = 0
    log10Pby4Refbase['T'] = 0
    for truebase in ['A', 'C', 'G', 'T']:
        prob = 0.0
        for pcrbase in ['A', 'C', 'G', 'T']:
            log10pcrbase = LOG10_1_MINUS_PCR_ERR
            if truebase != pcrbase:
                log10pcrbase = LOG10_PCR_ERR_PERSTATE
            if obsqual != 0:
                log10pcrbase = log10pcrbase + log10PofObsbaseGivenRefbase(obsbase, obsqual, pcrbase)
            prob = prob + math.pow(10, log10pcrbase)

        log10Pby4Refbase[truebase] = math.log10(prob)

    return log10Pby4Refbase


def calLLgenoByOneObs(obsbase, obsqual):
    log10Pby4Refbase = log10FourRefbase(obsbase, obsqual)
    LLgeno = initGenoLL()
    for gt in LLgeno:
        p_base = 0.0
        for g in gt:
            p_base = p_base + math.pow(10, log10Pby4Refbase[g.upper()]) / PLOIDY

        LLgeno[gt] = LLgeno[gt] + math.log10(p_base)

    return LLgeno


def calLLgenoDict(extbasemap, doRound=False):
    LLgenoDict = {}
    for site in extbasemap:
        if len(extbasemap[site]) > 0:
            LLgenoDict[site] = initGenoLL()
            for pair in extbasemap[site]:
                onellgeno = calLLgenoByOneObs(obsbase=pair[0], obsqual=pair[1])
                for each in onellgeno:
                    LLgenoDict[site][each] = LLgenoDict[site][each] + onellgeno[each]

    if doRound:
        for site in LLgenoDict:
            for each in LLgenoDict[site]:
                LLgenoDict[site][each] = round(LLgenoDict[site][each], 2)

    return LLgenoDict


class readPileupper():

    def __init__(self):
        self.basemap = {}
        self.fullbasemap = {}
        self.phasemap = {}
        self.extbasemap = {}

    def pileupBySample(self, readDict, SNPdb):
        if isinstance(readDict, dict):
            for mykey in readDict:
                self.pileupByRead(readDict[mykey], SNPdb)
                self.pileupFullBaseMapByRead(readDict[mykey])

        elif isinstance(readDict, list):
            for read in readDict:
                self.pileupByRead(read, SNPdb)
                self.pileupFullBaseMapByRead(read)

    def updateBySites(self, sites):
        rmsites = self.basemap.keys().difference(sites)
        for each in rmsites:
            del self.basemap[each]
            del self.phasemap[each]
            del self.extbasemap[each]
            del self.fullbasemap[each]

    def pileupFullBaseMapByRead(self, OneSAMrecordList):
        snpholder = {}
        sortedsnps = []
        read1 = OneSAMrecordList[0]
        start1 = read1.getAlignStart()
        stop1 = read1.getAlignEnd()
        start2 = start1
        stop2 = stop1
        snps1 = sets.Set(range(read1.getAlignStart(), read1.getAlignEnd() + 1))
        snpholder['1'] = snps1
        snpholder['2'] = sets.Set([])
        snpholder['ol'] = sets.Set([])
        sortedsnps = sorted(snps1)
        if len(OneSAMrecordList) > 1:
            read2 = OneSAMrecordList[1]
            start2 = read2.getAlignStart()
            stop2 = read2.getAlignEnd()
            snps2 = sets.Set(range(read2.getAlignStart(), read2.getAlignEnd() + 1))
            sortedsnps = sorted(snps1.union(snps2))
            olsnps = snps1.intersection(snps2)
            snpholder['ol'] = olsnps
            snpholder['2'] = snps2.difference(olsnps)
            snpholder['1'] = snps1.difference(olsnps)
        edges = [
         start1, stop1, start2, stop2]
        for i in xrange(len(sortedsnps)):
            sitei = sortedsnps[i]
            if sitei in snpholder['1']:
                c1 = read1.getCharAt(sitei)
                if isValidBase(c1):
                    if sitei not in self.fullbasemap:
                        self.fullbasemap[sitei] = 0
                    self.fullbasemap[sitei] = self.fullbasemap[sitei] + 1
            elif sitei in snpholder['2']:
                c1 = read2.getCharAt(sitei)
                if isValidBase(c1):
                    if sitei not in self.fullbasemap:
                        self.fullbasemap[sitei] = 0
                    self.fullbasemap[sitei] = self.fullbasemap[sitei] + 1
            else:
                c11 = read1.getCharAt(sitei)
                c12 = read2.getCharAt(sitei)
                useread = 0
                if read1.getAlignStart() == read2.getAlignStart() and len(read1.seq) != len(read2.seq):
                    if len(read2.getFormattedRead()) > len(read1.getFormattedRead()):
                        useread = 2
                    else:
                        useread = 1
                if isValidBase(c11) and (useread == 0 or useread == 1):
                    if sitei not in self.fullbasemap:
                        self.fullbasemap[sitei] = 0
                    self.fullbasemap[sitei] = self.fullbasemap[sitei] + 1
                if isValidBase(c12) and (useread == 0 or useread == 2):
                    if sitei not in self.fullbasemap:
                        self.fullbasemap[sitei] = 0
                    self.fullbasemap[sitei] = self.fullbasemap[sitei] + 1

    def pileupByRead(self, OneSAMrecordList, SNPdb):
        snpholder = {}
        sortedsnps = []
        read1 = OneSAMrecordList[0]
        start1 = read1.getAlignStart()
        stop1 = read1.getAlignEnd()
        start2 = start1
        stop2 = stop1
        snps1 = SNPdb.getSNPsitesNoDWithin(read1.getAlignStart(), read1.getAlignEnd()).intersection(SNPdb.withinRegionSites)
        snpholder['1'] = snps1
        snpholder['2'] = sets.Set([])
        snpholder['ol'] = sets.Set([])
        sortedsnps = sorted(snps1)
        if len(OneSAMrecordList) > 1:
            read2 = OneSAMrecordList[1]
            start2 = read2.getAlignStart()
            stop2 = read2.getAlignEnd()
            snps2 = SNPdb.getSNPsitesNoDWithin(read2.getAlignStart(), read2.getAlignEnd()).intersection(SNPdb.withinRegionSites)
            sortedsnps = sorted(snps1.union(snps2))
            olsnps = snps1.intersection(snps2)
            snpholder['ol'] = olsnps
            snpholder['2'] = snps2.difference(olsnps)
            snpholder['1'] = snps1.difference(olsnps)
        edges = [
         start1, stop1, start2, stop2]
        for i in xrange(len(sortedsnps)):
            sitei = sortedsnps[i]
            if sitei in snpholder['1']:
                c1 = read1.getCharAt(sitei)
                if isValidBase(c1):
                    if sitei not in self.basemap:
                        self.basemap[sitei] = initBaseMapBySite()
                    if sitei not in self.extbasemap:
                        self.extbasemap[sitei] = []
                    self.basemap[sitei][c1] = self.basemap[sitei][c1] + 1
                    self.basemap[sitei]['total'] = self.basemap[sitei]['total'] + 1
                    self.basemap[sitei][('to' + c1)].append(read1.qname)
                    self.extbasemap[sitei].append((c1, read1.getQualAt(sitei)))
            elif sitei in snpholder['2']:
                c1 = read2.getCharAt(sitei)
                if isValidBase(c1):
                    if sitei not in self.basemap:
                        self.basemap[sitei] = initBaseMapBySite()
                    if sitei not in self.extbasemap:
                        self.extbasemap[sitei] = []
                    self.basemap[sitei][c1] = self.basemap[sitei][c1] + 1
                    self.basemap[sitei]['total'] = self.basemap[sitei]['total'] + 1
                    self.basemap[sitei][('to' + c1)].append(read2.qname)
                    self.extbasemap[sitei].append((c1, read2.getQualAt(sitei)))
            else:
                c11 = read1.getCharAt(sitei)
                c12 = read2.getCharAt(sitei)
                useread = 0
                if read1.getAlignStart() == read2.getAlignStart() and len(read1.seq) != len(read2.seq):
                    if len(read2.getFormattedRead()) > len(read1.getFormattedRead()):
                        useread = 2
                    else:
                        useread = 1
                if isValidBase(c11) and (useread == 0 or useread == 1):
                    if sitei not in self.basemap:
                        self.basemap[sitei] = initBaseMapBySite()
                    if sitei not in self.extbasemap:
                        self.extbasemap[sitei] = []
                    self.basemap[sitei][c11] = self.basemap[sitei][c11] + 1
                    self.basemap[sitei]['total'] = self.basemap[sitei]['total'] + 1
                    self.basemap[sitei][('to' + c11)].append(read1.qname)
                    self.extbasemap[sitei].append((c11, read1.getQualAt(sitei)))
                if isValidBase(c12) and (useread == 0 or useread == 2):
                    if sitei not in self.basemap:
                        self.basemap[sitei] = initBaseMapBySite()
                    if sitei not in self.extbasemap:
                        self.extbasemap[sitei] = []
                    self.basemap[sitei][c12] = self.basemap[sitei][c12] + 1
                    self.basemap[sitei]['total'] = self.basemap[sitei]['total'] + 1
                    self.basemap[sitei][('to' + c12)].append(read2.qname)
                    self.extbasemap[sitei].append((c12, read2.getQualAt(sitei)))
                c1 = '0'
                if c11 == c12:
                    c1 = c11
            if isValidBase(c1) and sitei not in edges:
                for j in xrange(i + 1, len(sortedsnps)):
                    sitej = sortedsnps[j]
                    if sitej not in edges:
                        if sitej in snpholder['1']:
                            c2 = read1.getCharAt(sitej)
                        elif sitej in snpholder['2']:
                            c2 = read2.getCharAt(sitej)
                        else:
                            c2 = '0'
                            if read1.getCharAt(sitej) == read2.getCharAt(sitej):
                                c2 = read1.getCharAt(sitej)
                        if isValidBase(c2):
                            if (
                             sitei, sitej) not in self.phasemap:
                                self.phasemap[(sitei, sitej)] = initPhaseMapBySites()
                            self.phasemap[(sitei, sitej)][(c1, c2)] = self.phasemap[(sitei, sitej)][(c1, c2)] + 1
                            self.phasemap[(sitei, sitej)]['total'] = self.phasemap[(sitei, sitej)]['total'] + 1
                            break


class LLmaster():

    def __init__(self, LLgenoDict, phasemap, allelemaster, hladb, snpdb):
        self.LLstore = {}
        self.MaxLL = {}
        self.MaxLLpair = {}
        for locus in hladb.getLocus2Indices().keys():
            self.evaluateAllPairs(locus, LLgenoDict, phasemap, allelemaster, hladb, snpdb)

    def outputResults(self, outfile):
        if outfile is not None:
            filehandle = open(outfile, 'w')
            filehandle.write('Locus\tAllele1\tAllele2\tLLtot\tLLgeno\tLLphase\tLLfreq\n')
            for locus in sorted(self.MaxLL.keys()):
                for pair in self.MaxLLpair[locus]:
                    filehandle.write(locus + '\t' + pair[0] + '\t' + pair[1] + '\t' + str(self.MaxLL[locus]) + '\t' + str(self.LLstore[pair]['geno']) + '\t' + str(self.LLstore[pair]['phase']) + '\t' + str(self.LLstore[pair]['freq']) + '\n')

            filehandle.close()
        else:
            for locus in sorted(self.MaxLL.keys()):
                for pair in self.MaxLLpair[locus]:
                    print locus + '\t' + pair[0] + '\t' + pair[1] + '\t' + str(self.MaxLL[locus]) + '\t' + str(self.LLstore[pair]['geno']) + '\t' + str(self.LLstore[pair]['phase']) + '\t' + str(self.LLstore[pair]['freq'])

        return

    def evaluateAllPairs(self, locus, LLgenoDict, phasemap, allelemaster, hladb, snpdb):
        if len(locus.split('_')) < 2:
            locus = 'HLA_' + locus
        indices = sorted(hladb.Locus2Indices[locus])
        start = hladb.getStartPositions()[indices[0]]
        stop = hladb.getStopPositions()[indices[0]]
        positions = sorted(list(snpdb.getSNPsitesNoDWithin(start + 1, stop - 1).intersection(snpdb.withinRegionSites)))
        for i in xrange(len(indices)):
            aindex1 = indices[i]
            for j in xrange(i, len(indices)):
                aindex2 = indices[j]
                self.evaluateOnePair(aindex1, aindex2, positions, LLgenoDict, phasemap, allelemaster, hladb, snpdb)

    def evaluateOnePair(self, aindex1, aindex2, positionsForPhase, LLgenoDict, phasemap, allelemaster, hladb, snpdb):
        LLgeno = calAlleleLLgeno(aindex1, aindex2, LLgenoDict, hladb, snpdb)
        LLphase = calAlleleLLphase(aindex1, aindex2, positionsForPhase, phasemap, hladb, snpdb)
        aname1 = hladb.getNames()[aindex1]
        aname2 = hladb.getNames()[aindex2]
        LLfreq = math.log10(hladb.getMaxFreqs()[convertName(aname1, '4d')]) + math.log10(hladb.getMaxFreqs()[convertName(aname2, '4d')])
        pair = tuple(sorted([aname1, aname2]))
        pair4d = [ convertName(x, '4d') for x in list(pair) ]
        if LLgeno < 0:
            self.LLstore[pair] = {}
            self.LLstore[pair]['geno'] = LLgeno
            self.LLstore[pair]['phase'] = LLphase
            self.LLstore[pair]['freq'] = LLfreq
            self.LLstore[pair]['total'] = LLgeno + LLphase + LLfreq
            self.LLstore[pair]['pval'] = [allelemaster.pvals[pair4d[0]], allelemaster.pvals[pair4d[1]]]
            locus = convertName(aname1, 'locus')
            if locus not in self.MaxLL:
                self.MaxLL[locus] = self.LLstore[pair]['total']
            if locus not in self.MaxLLpair:
                self.MaxLLpair[locus] = []
            if self.LLstore[pair]['total'] > self.MaxLL[locus]:
                self.MaxLL[locus] = self.LLstore[pair]['total']
                self.MaxLLpair[locus] = [pair]
            elif self.LLstore[pair]['total'] == self.MaxLL[locus]:
                self.MaxLLpair[locus].append(pair)


def getResolution(pair, totlen=False):
    for rs in ['6d', '4d']:
        a1 = sets.Set([ convertName(x[0], rs) for x in pair ])
        if len(a1) == 1:
            break

    for rs in ['6d', '4d']:
        a2 = sets.Set([ convertName(x[1], rs) for x in pair ])
        if len(a2) == 1:
            break

    if totlen:
        return len(list(a1)[0].split(':')) + len(list(a2)[0].split(':'))
    else:
        if len(a1) == 1 and len(a2) == 1:
            return [list(a1)[0], list(a2)[0]]
        if len(a1) > 1:
            print a1
            return [
             (',').join(a1), list(a2)[0]]
        return [list(a1)[0], (',').join(a2)]


def getBestResolution(tuplelist):
    pair = [
     tuplelist[0]]
    for add in range(1, len(tuplelist)):
        tmp1 = copy.deepcopy(pair)
        tmp2 = copy.deepcopy(pair)
        tmp1.append(tuplelist[add])
        tmp2.append((tuplelist[add][1], tuplelist[add][0]))
        if getResolution(tmp2, True) > getResolution(tmp1, True):
            pair = tmp2
        else:
            pair = tmp1

    return getResolution(pair)


def outputResults(outfile, llmaster):
    filehandle = open(outfile, 'w')
    filehandle.write('Locus\tAllele1\tAllele2\tLLtot\tpval1\tpval2\n')
    for locus in sorted(llmaster.MaxLL.keys()):
        maxp1 = 0
        maxp2 = 0
        for pair in llmaster.MaxLLpair[locus]:
            p1 = llmaster.LLstore[pair]['pval'][0]
            p2 = llmaster.LLstore[pair]['pval'][1]
            if p1 == 'NA' and p2 == 'NA':
                maxp1 = 'NA'
                maxp2 = 'NA'
            else:
                if p1 > maxp1:
                    maxp1 = p1
                if p2 > maxp2:
                    maxp2 = p2

        if maxp1 == 'NA' and maxp2 == 'NA':
            filehandle.write(locus + '\t' + 'no call due to insufficient reads at this locus\n')
        else:
            if maxp1 < 1e-323:
                maxp1 = '0'
            elif isinstance(maxp1, float) and maxp1 > 0.01:
                maxp1 = '%.1e' % maxp1
            elif isinstance(maxp1, float):
                maxp1 = '%.1e' % maxp1
            if maxp2 < 1e-323:
                maxp2 = '0'
            elif isinstance(maxp2, float) and maxp2 > 0.01:
                maxp2 = '%.1e' % maxp2
            elif isinstance(maxp2, float):
                maxp2 = '%.1e' % maxp2
            if len(llmaster.MaxLLpair[locus]) > 1:
                pair = getBestResolution(llmaster.MaxLLpair[locus])
            else:
                pair = list(llmaster.MaxLLpair[locus][0])
            outpair0 = pair[0].split(',')
            if len(outpair0) == 1:
                outpair0 = convertName(pair[0], '', False)
            else:
                outpair0 = (',').join([ convertName(x, '', False) for x in outpair0 ])
            outpair1 = pair[1].split(',')
            if len(outpair1) == 1:
                outpair1 = convertName(pair[1], '', False)
            else:
                outpair1 = (',').join([ convertName(x, '', False) for x in outpair1 ])
            filehandle.write(locus + '\t' + outpair0 + '\t' + outpair1 + '\t' + str(round(llmaster.MaxLL[locus], 2)) + '\t' + str(maxp1) + '\t' + str(maxp2) + '\n')

    filehandle.close()


def calAlleleLLgeno(aindex1, aindex2, LLgenoDict, hladb, snpdb):
    allele1 = hladb.getHLAseqs()[aindex1]
    allele2 = hladb.getHLAseqs()[aindex2]
    start1 = hladb.getStartPositions()[aindex1]
    stop1 = hladb.getStopPositions()[aindex1]
    start2 = hladb.getStartPositions()[aindex2]
    stop2 = hladb.getStopPositions()[aindex2]
    llgeno = 0.0
    positions = snpdb.getSNPsitesNoDWithin(max(start1, start2) + 1, min(stop1, stop2) - 1).intersection(LLgenoDict.keys()).intersection(snpdb.withinRegionSites)
    for pos in positions:
        mykey = tuple(sorted([allele1[(pos - start1)], allele2[(pos - start2)]]))
        if mykey in LLgenoDict[pos]:
            llgeno = llgeno + LLgenoDict[pos][mykey]

    return llgeno


def calAlleleLLphase(aindex1, aindex2, positions, phasemap, hladb, snpdb):
    allele1 = hladb.getHLAseqs()[aindex1]
    allele2 = hladb.getHLAseqs()[aindex2]
    start1 = hladb.getStartPositions()[aindex1]
    stop1 = hladb.getStopPositions()[aindex1]
    start2 = hladb.getStartPositions()[aindex2]
    stop2 = hladb.getStopPositions()[aindex2]
    llphase = 0.0
    numInPhase = 0.0
    sumObs = 0.0
    sumInPhase = 0.0
    for i in xrange(len(positions) - 1):
        posi = positions[i]
        c11 = allele1[(posi - start1)]
        c21 = allele2[(posi - start2)]
        for j in xrange(i + 1, len(positions)):
            posj = positions[j]
            if (posi, posj) in phasemap:
                c12 = allele1[(posj - start1)]
                c22 = allele2[(posj - start2)]
                ishomo = True
                if c11 == c21 and c12 == c22:
                    numInPhase = phasemap[(posi, posj)][(c11, c12)]
                else:
                    ishomo = False
                    numInPhase = phasemap[(posi, posj)][(c11, c12)] + phasemap[(posi, posj)][(c21, c22)]
                numOutPhase = phasemap[(posi, posj)]['total'] - numInPhase
                addphaseLL = 0.0
                if ishomo:
                    addphaseLL = numInPhase * L_correct + numOutPhase * L_err
                else:
                    addphaseLL = numInPhase * math.log10((P_correct + P_err / 15.0) / 2.0) + numOutPhase * L_err
                llphase = llphase + addphaseLL
                break

    return llphase


def getGood(HLAhash, locus, qcut):
    out = sets.Set([])
    counter = 0
    tmphash = HLAhash[locus]
    length = len(tmphash)
    sortedv = sorted(tmphash.values())
    topQ = sortedv[int(len(sortedv) * qcut)]
    for allele4d in tmphash.keys():
        if tmphash[allele4d] >= topQ:
            out.add('HLA_' + allele4d)

    return out


def readMapFile(mapfile):
    HLAhash = {}
    filehandle = open(mapfile, 'r')
    for strline in filehandle:
        strline = strline.strip()
        parts = strline.split('*')
        locus = parts[0]
        allele = strline
        allele = convertName(allele, '4d')
        if locus not in HLAhash:
            HLAhash[locus] = {}
            HLAhash[locus][allele] = 1
        elif allele not in HLAhash[locus]:
            HLAhash[locus][allele] = 1
        else:
            HLAhash[locus][allele] = HLAhash[locus][allele] + 1

    filehandle.close()
    return HLAhash


def findAllelesToSearch(mapfile=None):
    mapper4d = {}
    HLAhash = {}
    allelesToSearch = sets.Set([])
    if mapfile is None:
        return allelesToSearch
    else:
        filehandle = open(mapfile, 'r')
        for strline in filehandle:
            strline = strline.strip()
            parts = strline.split('*')
            locus = parts[0]
            allele = strline
            allele4d = convertName(allele, '4d')
            if allele4d not in mapper4d:
                mapper4d[allele4d] = sets.Set([])
            mapper4d[allele4d].add(allele)
            if locus not in HLAhash:
                HLAhash[locus] = {}
                HLAhash[locus][allele] = 1
            elif allele not in HLAhash[locus]:
                HLAhash[locus][allele] = 1
            else:
                HLAhash[locus][allele] = HLAhash[locus][allele] + 1

        filehandle.close()
        loci = sorted(HLAhash.keys())
        for locus in loci:
            counter = 0
            tmphash = HLAhash[locus]
            length = len(tmphash)
            sortedv = sorted(tmphash.values())
            quantile = 0.7
            if locus == 'A' or locus == 'B' or locus == 'C' or locus == 'DRB1':
                quantile = 0.9
            topQ = sortedv[int(len(sortedv) * quantile)]
            if topQ <= 5 and len(sortedv) <= 100:
                topQ = 1
            for allele in tmphash.keys():
                if tmphash[allele] >= topQ:
                    allele4d = convertName(allele, '4d')
                    allelesToSearch = allelesToSearch.union(mapper4d[allele4d])

        return allelesToSearch


def getRAllele1(alist, count0, hladb, maptable):
    out = sets.Set([])
    afreq = [ hladb.getMaxFreqs()[x] for x in alist ]
    if len(alist) > 1 and min(afreq) < max(afreq) and count0 < 400:
        locus = convertName(alist[0], resolution='locus', verbal=False)
        goodones = getGood(maptable, locus, 0.95)
        for x in alist:
            if x not in goodones:
                out.add(x)

    return out


def getRAllele0(alist, count0, hladb, maptable):
    out = sets.Set([])
    afreq = [ hladb.getMaxFreqs()[x] for x in alist ]
    if len(alist) > 1 and min(afreq) < max(afreq) and count0 < 900:
        locus = convertName(alist[0], resolution='locus', verbal=False)
        goodones = getGood(maptable, locus, 0.97)
        for x in alist:
            if x not in goodones:
                out.add(x)

    return out


def initBaseMapBySite():
    check = {}
    check['A'] = 0
    check['T'] = 0
    check['C'] = 0
    check['G'] = 0
    check['total'] = 0
    check['toA'] = []
    check['toC'] = []
    check['toG'] = []
    check['toT'] = []
    return check


def initPhaseMapBySites():
    check = {}
    check[('A', 'A')] = 0
    check[('C', 'C')] = 0
    check[('G', 'G')] = 0
    check[('T', 'T')] = 0
    check[('A', 'C')] = 0
    check[('A', 'G')] = 0
    check[('A', 'T')] = 0
    check[('C', 'G')] = 0
    check[('C', 'T')] = 0
    check[('G', 'T')] = 0
    check[('C', 'A')] = 0
    check[('G', 'A')] = 0
    check[('T', 'A')] = 0
    check[('G', 'C')] = 0
    check[('T', 'C')] = 0
    check[('T', 'G')] = 0
    check['total'] = 0
    return check


def initGenoLL():
    LLgeno = {}
    LLgeno[('A', 'C')] = 0.0
    LLgeno[('A', 'G')] = 0.0
    LLgeno[('A', 'T')] = 0.0
    LLgeno[('C', 'G')] = 0.0
    LLgeno[('C', 'T')] = 0.0
    LLgeno[('G', 'T')] = 0.0
    LLgeno[('A', 'A')] = 0.0
    LLgeno[('C', 'C')] = 0.0
    LLgeno[('G', 'G')] = 0.0
    LLgeno[('T', 'T')] = 0.0
    return LLgeno


def initCoverRegion(mytype='CDS'):
    regions = []
    if mytype.lower() == 'cds':
        regions.append([29910331, 29910403])
        regions.append([29910534, 29910803])
        regions.append([29911045, 29911320])
        regions.append([29911899, 29912174])
        regions.append([29912277, 29912393])
        regions.append([29912836, 29912868])
        regions.append([29913011, 29913058])
        regions.append([29913228, 29913232])
        regions.append([31322260, 31322303])
        regions.append([31322410, 31322442])
        regions.append([31322884, 31323000])
        regions.append([31323094, 31323369])
        regions.append([31323944, 31324219])
        regions.append([31324465, 31324734])
        regions.append([31324863, 31324935])
        regions.append([31239776, 31239848])
        regions.append([31237115, 31237162])
        regions.append([31237270, 31237302])
        regions.append([31237743, 31237862])
        regions.append([31237987, 31238262])
        regions.append([31238850, 31239125])
        regions.append([31239376, 31239645])
        regions.append([31236946, 31236950])
        regions.append([32605236, 32605317])
        regions.append([32609087, 32609335])
        regions.append([32609749, 32610030])
        regions.append([32610387, 32610541])
        regions.append([32628013, 32628026])
        regions.append([32628636, 32628659])
        regions.append([32629124, 32629234])
        regions.append([32629744, 32630025])
        regions.append([32632575, 32632844])
        regions.append([32634276, 32634384])
        regions.append([32546868, 32546881])
        regions.append([32548024, 32548047])
        regions.append([32548523, 32548633])
        regions.append([32549334, 32549615])
        regions.append([32551886, 32552155])
        regions.append([32557420, 32557519])
    elif mytype.lower() == 'exon':
        regions.append([29910247, 29910403])
        regions.append([29910534, 29910803])
        regions.append([29911045, 29911320])
        regions.append([29911899, 29912174])
        regions.append([29912277, 29912393])
        regions.append([29912836, 29912868])
        regions.append([29913011, 29913058])
        regions.append([29913228, 29913532])
        regions.append([31322256, 31322303])
        regions.append([31322410, 31322442])
        regions.append([31322884, 31323000])
        regions.append([31323094, 31323369])
        regions.append([31323944, 31324219])
        regions.append([31324465, 31324734])
        regions.append([31324863, 31324989])
        regions.append([31236775, 31239950])
        regions.append([31237115, 31237162])
        regions.append([31237270, 31237302])
        regions.append([31237743, 31237862])
        regions.append([31237987, 31238262])
        regions.append([31238850, 31239125])
        regions.append([31239376, 31239645])
        regions.append([31239776, 31239913])
        regions.append([32605183, 32605317])
        regions.append([32609087, 32609335])
        regions.append([32609749, 32610030])
        regions.append([32610387, 32610561])
        regions.append([32610541, 32610987])
        regions.append([32627808, 32628026])
        regions.append([32628636, 32628659])
        regions.append([32629124, 32629234])
        regions.append([32629744, 32630025])
        regions.append([32632575, 32632844])
        regions.append([32634276, 32634466])
        regions.append([32546547, 32546881])
        regions.append([32548024, 32548047])
        regions.append([32548523, 32548633])
        regions.append([32549334, 32549615])
        regions.append([32551886, 32552155])
        regions.append([32557420, 32557613])
    else:
        print 'ERR: please use either CDS or EXON as coverage region'
    return regions


def countBaseByRegions(regions):
    bcounts = sum([ x[1] - x[0] + 1 for x in regions ])
    return bcounts


def getBaseByRegions(regions):
    sites = []
    for x in regions:
        sites.extend(range(x[0], x[1] + 1))

    return sites


def estimateCov(basemap, method, addmissed):
    sites = getBaseByRegions(initCoverRegion('CDS'))
    ol = sets.Set(sites).intersection(sets.Set(basemap.keys()))
    missed = sets.Set(sites).difference(sets.Set(basemap.keys()))
    kept = [ basemap[x] for x in ol ]
    if len(missed) > 0 and addmissed:
        kept.extend([ 0 for x in missed ])
    kept = sorted(kept)
    if method == 'median':
        if not len(kept) % 2:
            return (kept[(len(kept) / 2)] + kept[(len(kept) / 2 - 1)]) / 2.0
        return kept[(len(kept) / 2)]
    else:
        return sum(kept) * 1.0 / len(kept)