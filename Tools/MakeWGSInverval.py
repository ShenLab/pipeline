#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# MakeWGSInverval.py
#=========================================================================

import argparse
import csv
import re

contig_to_split = [str(i) for i in xrange(1, 23)]
contig_to_split.extend(['X','Y'])
print "congig to split:",contig_to_split
contig_to_split = set(contig_to_split)

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fai', type=str, help='fai file contains genomes contig and length')
    parser.add_argument('-l', '--length', type=int, default=10000, help='target interval length. Default 10000')
    parser.add_argument('-o', '--out', type=str, default="WGS.bed", help='OutName for Interval bed file')
    args = parser.parse_args()
    return args.fai, args.length, args.out

class MakeWGSInverval:
    def __init__(self, fai, length, outname):
        self.fai = fai
        self.length = length
        self.outname = outname
    def run(self):
        fout = open(self.outname, 'wb') 
        fout.write("#Contig\tStart\tEnd\n")
        contigs = csv.reader(open(self.fai, 'rb'), delimiter="\t")
        for row in contigs:
            print row
            contig_name, contig_length = row[0], row[1]
            if contig_name in contig_to_split:
                for start, end in self.split(contig_name, contig_length):
                    fout.write("{}\t{}\t{}\n".format(str(contig_name), str(start), str(end)))
            else:
                fout.write("{}\t{}\t{}\n".format(str(contig_name), str(0), str(contig_length)))
        fout.close()
    def split(self, contig_name, contig_length):
        Start = 0
        End = self.length
        while int(End) < int(contig_length):
            yield Start, End
            Start = End + 1
            End = min(End + self.length, int(contig_length))
        yield Start, End

def main():
    fai, length, outname = GetOptions()
    ins = MakeWGSInverval(fai, length, outname)
    ins.run()
    return


if __name__ == '__main__':
    main()
