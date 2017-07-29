#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# ExmAln.0a.Confirm_SampleID_from_Fastq.py
#=========================================================================

import argparse
import re

#=========================================================================
# Modify this part each time run. Customized by different Batch of Data
SM = '(OMG\d+-\d+-[A-Za-z0-9-]+)'
ID = '_([A-Za-z0-9-]+)_'
PL = 'Illumina'
LIB = 'LIB'
PAIR = '_R(1|2)'
CN = 'OMG'
#=========================================================================

name = re.compile(SM) #PIPseq ID
rgid = re.compile(ID)

class Sample:
    def __init__(self, ID, F_name, F_path):
        self.SampleID = ID
        #self.RGID = [rgid.search(F_name).group(1)]
        self.RGID = [F_path]
        self.F_paths = [F_path]
    def addPath(self, F_path):
        self.F_paths.append(F_path)
        F_name = F_path.split('/')[-1]
        #self.RGID.append(rgid.search(F_name).group(1))
        self.RGID.append(F_path)
        #print F_name, rgid.search(F_name).group(1)

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq', type=str, help='Fastq List')
    args = parser.parse_args()

    return args.fastq

def ExtractSample(FastqList):
    res = {}
    fin = open(FastqList, 'rb')
    for l in fin:
        fpath = l.strip()
        fname = fpath.split('/')[-1]
        sample = name.search(fname).group(1)
        #print sample
        if sample in res:
            res[sample].addPath(fpath)
        else:
            res[sample] = Sample(sample, fname, fpath)
    fin.close()
    for sample in sorted(res.values(),key=lambda x:x.SampleID):
        print sample.SampleID
        print '\n'.join(sample.RGID)

def MakeNextStep(FastqList):
    StepName = 'run_MakeFastqTable.sh'
    fout = open(StepName, 'wb')
    fout.write('#!/bin/bash\n')
    fout.write('InpFil={}\n'.format(FastqList))
    fout.write('CMD={}\n'.format('/home/yufengshen/CUMC/Exome-pipeline-Jiayao/ExmAln.0b.MakeFastqTable.py'))
    fout.write('Target=\n'.format())
    fout.write('\n'.format())

def main():
    FastqList = GetOptions()
    ExtractSample(FastqList)
    #MakeNextStep()
    return


if __name__ == '__main__':
    main()
