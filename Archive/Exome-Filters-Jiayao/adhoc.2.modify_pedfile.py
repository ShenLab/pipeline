#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Modify Ped files
#=========================================================================

from optparse import OptionParser


def GetOptions():
    parser = OptionParser()
    parser.add_option('-p', '--ped', dest='PEDFILE',
                      metavar='PEDFILE', help='Error Pedigree File to be Modify')
    parser.add_option('-o', '--out', dest='OUTFILE', metavar='OUTFILE',
                      help='Output File Name of Modified Ped File')
    (options, args) = parser.parse_args()
    if options.OUTFILE == None:
        options.OUTFILE = '.'.join(options.PEDFILE.split('.')[
                                   :-1]) + '_Modifed.ped'
    return options.PEDFILE, options.OUTFILE


def ModifyPed(pedfile, outname):
    fin = open(pedfile, 'rb')
    fout = open(outname, 'wb')
    fin.readline()
    fout.write('#FamID\tProband\tFatherID\tMotherID\tGender\tPhenotype\n')
    for l in fin:
        Proband, father, mother, phenopb, phenoFA, phenoMA = l.strip().split('\t')
        fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            Proband, Proband, father, mother, '0', phenopb))
        fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            Proband, father, '0', '0', '1', phenoFA))
        fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            Proband, mother, '0', '0', '2', phenoMA))


def main():
    pedfile, outname = GetOptions()
    ModifyPed(pedfile, outname)
    return


if __name__ == '__main__':
    main()
