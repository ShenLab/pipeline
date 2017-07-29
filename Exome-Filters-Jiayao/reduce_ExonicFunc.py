#!/home/local/users/jw/bin/python2.7
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
#
#=========================================================================

from optparse import OptionParser
import gzip
from sys import argv
import re


def GetOptions():
    parser = OptionParser()
    parser.add_option('-', '--', dest='', metavar='', help='')
    (options, args) = parser.parse_args()

    return


def main():
    if argv[1].endswith('.gz'):
        fin = gzip.open(argv[1], 'rb')
    else:
        fin = open(argv[1], 'rb')
    res_1 = {}
    res_2 = {}
    for l in fin:
        if not l.startswith('##'):
            llist = l.split('\t')
            INFO = llist[7]
            tmp_1 = re.findall('ExonicFunc.refGene=(.+?);', INFO)
            for item in tmp_1:
                if item in res_1:
                    res_1[item] += 1
                else:
                    res_1[item] = 1
            tmp_2 = re.findall(';Func.refGene=(.+?);', INFO)
            for item in tmp_2:
                if item in res_2:
                    res_2[item] += 1
                else:
                    res_2[item] = 1
    print res_1
    print res_2


if __name__ == '__main__':
    main()
