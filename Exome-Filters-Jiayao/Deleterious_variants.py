#!/home/local/users/jw/bin/python2.7
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Find deleterious variants in a vcf file
# Rare: ExAC_ALL < 1% (0.01)
# LGD: 'splicing', 'frameshift insertion','frameshift deletion','frameshift block substitution','stopgain','stoploss'
# D-mis: MetaSVM_pred=D
#=========================================================================

from optparse import OptionParser
from Variant import INFO
import gzip


def LGD_Dmis(info):
    for annotation in info.annovar:
        flag = False
        try:
            if float(annotation['ExAC_ALL']) >= 0.01:
                continue
        except ValueError:
            pass
        if annotation['ExonicFunc.refGene'] == 'nonsynonymous_SNV':
            if 'D' == annotation['MetaSVM_pred']:
                flag = True
            try:
                if float(annotation['CADD_phred']) > 20:  # D-mis
                    flag = True
            except ValueError:
                continue
        if ('splicing' in annotation['Func.refGene']) or (annotation['Func.refGene'] == 'exonic' and annotation['ExonicFunc.refGene'] in ['stoploss', 'stopgain', 'frameshift_insertion', 'frameshift_deletion', 'frameshift block substitution']):  # LGD
            flag = True
        if flag == True:
            return True
    return False


def GetOptions():
    parser = OptionParser()
    parser.add_option('-v', '--vcf', dest='VCF',
                      metavar='VCF', help='Input VCF')
    parser.add_option('-o', '--out', dest='OUT',
                      metavar='OUT', help='Output VCF')
    (options, args) = parser.parse_args()

    return options.VCF, options.OUT


def Filter_VCF(In, Out):
    if In.endswith('.gz'):
        fin = gzip.open(In, 'rb')
    else:
        fin = open(In, 'rb')
    if Out.endswith('.gz'):
        fout = gzip.open(Out, 'wb')
    else:
        fout = open(Out, 'wb')
    counter = 0
    for l in fin:
        if l.startswith('#'):
            fout.write(l)
        else:
            counter += 1
            llist = l.split('\t')
            info = INFO(llist[4].split(','), llist[7])
            if LGD_Dmis(info):
                fout.write(l)
        if counter % 100 == 0 and counter != 0:
            print 'processed %d Record.' % counter
    fin.close()
    fout.close()


def main():
    In, Out = GetOptions()
    Filter_VCF(In, Out)


if __name__ == '__main__':
    main()
