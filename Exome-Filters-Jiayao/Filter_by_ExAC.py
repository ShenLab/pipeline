#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Filters on VCF file
#=========================================================================

from optparse import OptionParser
import gzip
import os
import re

VQSR = re.compile('([\d.]+)')


def GetOptions():
    parser = OptionParser()
    parser.add_option("-v", "--vcf", dest="VCF",
                      help="Input VCF file", metavar="VCFfile")
    parser.add_option("-o", "--outvcf", dest="OutVCF",
                      help="Name of Output VCF file", metavar="OutVCF", default="Filterd.vcf")
    parser.add_option("-f", "--filters", dest="Filters", metavar="Filters",
                      help="Filters apply on the variants. Splited by ','"
                      )
    (options, args) = parser.parse_args()
    VCFin = options.VCF
    VCFout = options.OutVCF
    Filters = ""
    if options.Filters != None:
        Filters = options.Filters.split(',')

    return VCFin, VCFout, Filters


def Filter(VCFin, VCFout, Filters):
    if VCFin.endswith('.gz'):
        hand = gzip.open(VCFin, 'rb')
    else:
        hand = open(VCFin, 'rb')
    if VCFout.endswith('.gz'):
        fout = gzip.open(VCFout, 'wb')
    else:
        fout = open(VCFout, 'wb')
    Count_All = 0
    Count_Pass = 0
    for l in hand:
        if l.startswith('##'):
            fout.write(l)
            continue
        elif l.startswith('#'):
			fout.write('##FiltersByRareCoing. ExAC_ALL:1e-2. RefGene.Func:exonic,splicing,exonic-splicing\n')
			fout.write(l)
			continue
        llist = l.strip().split('\t')
        Count_All += 1

        #CodingRegion = ['exonic','splicing']
        CodingRegion = ['exonic', 'splicing', 'exonic-splicing']
        if F_ExAC_Coding(llist[7], 1e-2, CodingRegion):
            fout.write(l)
            Count_Pass += 1

        if Count_All % 10000 == 0:
            print "Read %d Variants" % Count_All
    fout.close()
    print "Finish Reading All %d Variants. %d Variants Pass the Filter" % (Count_All, Count_Pass)

def F_ExAC_Coding(INFO, ExACfreq, Region):
    infolist = INFO.split(';')
    infodict = {}
    for kv in infolist:
        kv = kv.split('=')
        if len(kv) == 2:
            k, v = kv
            if k in infodict:
                infodict[k] = infodict[k] + ',' + v
            else:
                infodict[k] = v

    Pass_region = False
    for item, v in zip(Region, infodict['Func.refGene'].split(',')):
        if item in v:
            Pass_region = True
            continue
    if Pass_region == False:
        return False

    for v in infodict['ExAC_ALL'].split(','):
        try:
            if float(v) <= ExACfreq:
                return True
        except ValueError:
            return True
    return False


def main():
    VCFin, VCFout, Filters = GetOptions()
    Filter(VCFin, VCFout, Filters)


if __name__ == '__main__':
    main()
