#!/home/local/users/jw/bin/python2.7
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Split a multi-fam ped file to different pedigrees
#=========================================================================

from optparse import OptionParser


def GetOptions():
    parser = OptionParser()
    parser.add_option('-p', '--ped', dest='PED',
                      metavar='PED', help='Pedigree file to split')
    (options, args) = parser.parse_args()

    return options.PED


def main():
    ped = GetOptions()
    fin = open(ped, 'rb')
    for l in fin:
        if l.startswith('#'):  # header line
            header = l
            last_fam = None
        else:
            llist = l.strip().split('\t')
            fam = llist[0]
            if fam != last_fam:
                if 'fout' in locals():
                    fout.close()
                last_fam = fam
                fout = open(fam + '.ped', 'wb')
                fout.write(header)
                fout.write(l)
            else:
                fout.write(l)
    fout.close()


if __name__ == '__main__':
    main()
