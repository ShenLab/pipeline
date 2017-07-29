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
    parser.add_option("-l", "--lowerBound", dest="lowerBound",
                      help="lowerBound AF for filtering")
    parser.add_option("-u", "--upperBound", dest="upperBound",
                      help="upperBound AF for filtering")
    (options, args) = parser.parse_args()
    if options.lowerBound == None:
        options.lowerBound == '0'
    if options.upperBound == None:
        options.upperBound == '1'
    return options.VCF, options.OutVCF, options.lowerBound, options.upperBound


def Filter(VCFin, VCFout, lower, upper):
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
            fout.write('##FiltersBy gnomAD freq {} - {}\n'.format(lower, upper))
            fout.write(l)
            continue
        llist = l.strip().split('\t')
        Count_All += 1
        lower = float(lower)
        upper = float(upper)
        if F_WGSAF(llist[7], lower, upper):
            fout.write(l)
            Count_Pass +=1
        

        if Count_All % 10000 == 0:
            print "Read %d Variants" % Count_All
    fout.close()
    print "Finish Reading All %d Variants. %d Variants Pass the Filter" % (Count_All, Count_Pass)


def F_WGSAF(INFO, lower, upper):
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
    AF = infodict['gnomAD_genome_ALL'].split(',')[-1]
    if AF == '.':
        AF = 0
    if lower == upper:
        if float(AF) == lower:
            return True
        else:
            return False
    else:
        if float(AF)  > lower and float(AF) <= upper:
            return True
        else:
            return False

    #for AF in infodict['gnomAD_genome_ALL'].split(','):
    #    try:
    #        if float(AF) <= cutoff:
    #            return True
    #    except ValueError:
    #        if AF == '.':
    #            return True
    #        else:
    #            print AF
    #            return False
    #return False


def main():
    VCFin, VCFout, lower, upper = GetOptions()
    Filter(VCFin, VCFout, lower, upper)


if __name__ == '__main__':
    main()
