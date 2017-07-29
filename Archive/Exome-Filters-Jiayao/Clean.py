#!/home/local/users/jw/anaconda2/bin/python
# -*- coding: utf-8 -*-
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
#
#=========================================================================
import gzip
from optparse import OptionParser


def has_var(samples):
    for sample in samples:
        GT = sample.split(':')[0]
        if GT != './.' and GT != '0/0':
            return True
    return False


def See(Count, num):
    if Count % num == 0:
        print "Read %d variants." % Count


def Clean(InpFil, OutFil):
    if InpFil.endswith(".gz"):
        fin = gzip.open(InpFil, 'rb')
    else:
        fin = open(InpFil, 'rb')
    if OutFil.endswith('.gz'):
        fout = gzip.open(OutFil, 'wb')
    else:
        fout = open(OutFil, 'wb')
    Count, rm = 0, 0
    for l in fin:
        if l.startswith("#"):
            fout.write(l)
            continue
        Count += 1
        samples = l.strip().split('\t')[9:]
        if has_var(samples):
            fout.write(l)
        else:
            rm += 1
        See(Count, 10000)

    print "Read %d Total, %d pass, %d removed." % (Count, Count - rm, rm)


def GetOptions():
    parser = OptionParser()
    parser.add_option('-v', '--vcf', dest='VCF', metavar='VCF',
                      help='Input VCF file to be clean')
    parser.add_option('-o', '--out', dest='OUT',
                      metavar='OUT', help='Output VCF file')
    (options, args) = parser.parse_args()
    return options.VCF, options.OUT


def main():
    InpFil, OutFil = GetOptions()
    Clean(InpFil, OutFil)


if __name__ == '__main__':
    main()
