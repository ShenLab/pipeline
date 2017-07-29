#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# ExmAdHoc.2a.MergeBams.py
#=========================================================================

import argparse
import re
import os 

SampleID = re.compile('(OMG[A-Za-z0-9-]+)') # Need to be Modified each time according to BamName convention.
REF = '/home/yufengshen/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b38.biocluster.sh'
CMD = '/home/yufengshen/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.2.MergeBams.sh'
Metrix = 'false'

class BAM:
    def __init__(self, FilPath):
        self.FilPath = FilPath.strip()
        self.FilName = self.FilPath.split('/')[-1]
        self.SampleID = SampleID.search(self.FilName).group(0)
        print self.SampleID, self.FilName

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamlist', type=str, help='Bam list contains a list of bam files.')
    parser.add_argument('-o', '--outname', type=str, help='Output File Name.')
    args = parser.parse_args()
    if args.outname == None:
        args.outname = 'run_'+args.bamlist.split('/')[-1].rstrip('.list').rstrip('.txt')+'.MergeBam.sh'
        #args.outname = re.search('',args.bamlist)
    return args.bamlist, args.outname

def ManipulateBamList(bamlist, outname):
    samples = {}
    fin = open(bamlist, 'rb')
    for l in fin:
        bam = BAM(l)
        #print bam.FilPath
        if bam.SampleID not in samples:
            samples[bam.SampleID] = [bam]
        else:
            #print bam.SampleID, samples[bam.SampleID]
            samples[bam.SampleID].append(bam)
    
    fout = open(outname,'wb')
    fout.write('#!/bin/bash\n')
    fout.write('REF={}\n'.format(REF))
    fout.write('CMD={}\n\n'.format(CMD))
    cwd = os.getcwd()
    for sample,bams in samples.items():
        print sample, bams
        if len(bams) == 1:
            fout.write('mv {}.bam {}.bam \n'.format(bams[0].FilPath.rstrip('.bam'), cwd+'/'+bams[0].SampleID))
            fout.write('mv {}.bai {}.bai \n'.format(bams[0].FilPath.rstrip('.bam'), cwd+'/'+bams[0].SampleID))
        else:
            TmpInput = '{}.bam.list'.format(bams[0].SampleID)
            with open(TmpInput,'wb') as tmpfout:
                for bam in bams:
                    tmpfout.write(bam.FilPath+'\n')
            fout.write('qsub $CMD -r $REF -i {} -m {}\n'.format(TmpInput, Metrix))
    fin.close()
    fout.close()
def main():
    bamlist, outname = GetOptions()
    ManipulateBamList(bamlist, outname)
    return


if __name__ == '__main__':
    main()
