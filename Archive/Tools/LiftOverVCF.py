#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# LiftOverVCF.py
#=========================================================================

import argparse
import subprocess
import gzip

Hg19ToHg38 = '/home/yufengshen/resources/reference_genomes/LiftOver/hg19ToHg38.over.chain'
Hg38ToHg19 = '/home/yufengshen/resources/reference_genomes/LiftOver/hg38ToHg19.over.chain'

class LiftOver:
    def __init__(self, Vcf, Out, Chain):
        #self.Chain = Hg19ToHg38
        self.Chain = Hg38ToHg19 
        self.OriginalVCF = Vcf
        self.OriginalHand = openVCF(self.OriginalVCF)
        self.OriginalBed = self.OriginalVCF + '.bed'
        self.OutVCF = Out
        self.OutBed = Out+ '.bed'
        self.OutUmappedBed = Out + '.Umapped.bed'
        self.OutHand = open(self.OutVCF, 'wb')
    def Vcf2Bed(self):
        print "converting original vcf into bed file"
        tmp_fout = open(self.OriginalBed,'wb')
        for l in self.OriginalHand:
            if l.startswith('#'):
                continue
            record = VCFRecord(l)
            tmp_fout.write(record.FormOut())
        print 'Done'
    def RunLiftOver(self):
        print '\nRunning LiftOver...'
        cmd = 'liftOver {} {} {} {}'.format(self.OriginalBed, self.Chain, self.OutBed, self.OutUmappedBed)
        liftover = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        liftover.wait()
        print 'Done'
        return liftover.returncode
    
    def Bed2Vcf(self):
        print '\nConvering liftovered bed into VCF'
        OutBed = open(self.OutBed, 'rb')
        self.OriginalHand.seek(0)
        for l in self.OriginalHand:
            if l.startswith('##'):
                self.OutHand.write(l)
            elif l.startswith('#'):
                self.OutHand.write('##PROGRAM<LiftOver From HG19 to HG38>\n')
                self.OutHand.write(l)
            else:
                break
        self.OriginalHand.close()
        for l in OutBed:
            record = BedRecord(l)
            self.OutHand.write(record.FormOut())
        self.OutHand.close()
        print 'Done'

class VCFRecord:
    def __init__(self, l):
        llist = l.strip().split('\t')
        self.Chr, self.Pos = llist[0], llist[1]
        self.Other = '%%'.join(llist[2:])
        self.Start = str(int(self.Pos) - 1)
        self.End = str(self.Pos)
        if not self.Chr.startswith('chr'):
            self.Chr = 'chr' + self.Chr
    def FormOut(self):
        return '{}\t{}\t{}\t{}\n'.format(self.Chr, self.Start, self.End, self.Other)

class BedRecord:
    def __init__(self, l):
        llist = l.strip().split('\t')
        self.Chr, self.Start, self.End = llist[0], llist[1], llist[2]
        self.Pos = self.End
        self.Other = llist[3].split('%%')
        if self.Chr.startswith('chr'):
            self.Chr = self.Chr[3:]
    def FormOut(self):
        return '{}\t{}\t{}\n'.format(self.Chr, self.Pos, '\t'.join(self.Other))

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-V', '--InputVCF', type=str, help='VCF file to be liftover')
    parser.add_argument('-N', '--NewVCF', type=str, help='Name of New VCF file after LiftOver')
    parser.add_argument('-C', '--Chain', type=str, help='Chain used for LiftOver')
    args = parser.parse_args()
    
    return args.InputVCF, args.NewVCF, args.Chain

def openVCF(VCFname):
    if VCFname.endswith('.vcf.gz'):
        return gzip.open(VCFname, 'rb')
    elif VCFname.endswith('.vcf'):
        return open(VCFname, 'rb')
    else:
        return open(VCFname, 'rb')
def main():
    Vcf, Out, Chain = GetOptions()
    lift =  LiftOver(Vcf, Out, Chain)
    lift.Vcf2Bed()
    lift.RunLiftOver()
    lift.Bed2Vcf()
    return


if __name__ == '__main__':
    main()
