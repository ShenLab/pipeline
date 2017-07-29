#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# PSAP_PrepareInput.py
# Prepare bamout script and variant table for IGV shot
# Input: Dir contain subdirs of PSAP results, (PSAP.csv, .ped file is need)
#        TXT file contains (sampleID,BAM) tuple
# Output: List of (sampleID Interval1 Interval2 ... Inverbaln)
#         Scripts to run BAMOUT.
#         Scripts to run IGV
#=========================================================================

import argparse
import csv
import re
import os

# SampleName = re.compile('(COL-CHUNG_[a-zA-z0-9-]+_Diabetes_[\d]+)*.bam')
# #Example for Sample name in RGN data
SampleName = re.compile('(CARE[a-zA-Z0-9-]+).bam')
ProjectName = re.compile('\w+')
EXM_REF = '$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.biocluster.sh'
BAMOUT_CMD = '$HOME/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.15b.BamOut.sh'
IGV_GENERATE_CMD = '$HOME/CUMC/Exome-IGV-Jiayao/Generate_IGV_plots.sh'


def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', type=str, required=True,
                        help='Directory contains PSAP by fam results, each subdir should have .csv and .ped')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='File contains bam locations. ')
    parser.add_argument('-l1', '--list1', type=str,
                        help='list of sample and corresponding intervals, For Bamout')
    parser.add_argument('-l2', '--list2', type=str,
                        help='list of variants and corresponding bam, For IGV')
    parser.add_argument('-s1', '--script1', type=str,
                        help='Script to run GATK Bamout. ')
    parser.add_argument('-s2', '--script2', type=str,
                        help='Script to run IGV. ')
    parser.add_argument('-p', '--parallel', type=int,
                        default=20, help='Number of Parallel to go. ')
    parser.add_argument('--popscore', type=float, default=1e-3,
                        help='popscore cutoff to generate. ')
    parser.add_argument('--chet', type=bool, default=True,
                        help='Whether to generate Chet candidates. ')
    parser.add_argument('--snp', type=bool, default=False,
                        help='Whether to generate snp candidates ')
    args = parser.parse_args()
    args.dir = os.path.abspath(args.dir)
    print args.dir
    if args.list1 == None:
        args.list1 = ProjectName.search(
            args.dir.split('/')[-1]).group(0) + '.samplelist'
    if args.list2 == None:
        args.list2 = ProjectName.search(
            args.dir.split('/')[-1]).group(0) + '.varlist'
    if args.script1 == None:
        args.script1 = 'run_bamout.sh'
    if args.script2 == None:
        args.script2 = 'run_IGV.sh'
    return args


class PSAP_PrepareIGVInput:
    # def __init__(self, Dir, BAMList, VarList1, VarList2, Script1, Script2):
    def __init__(self, args):
        self.Dir = os.path.abspath(args.dir)
        self.bampath = BamLocation(args.bam)
        self.VarList1 = os.path.abspath(args.list1)
        self.VarList1_hand = open(self.VarList1, 'wb')
        self.VarList2 = os.path.abspath(args.list2)
        self.VarList2_hand = open(self.VarList2, 'wb')
        self.Script1 = open(args.script1, 'wb')
        self.Script2 = open(args.script2, 'wb')
        self.popScore = args.popscore
        self.parallel = args.parallel
        self.chet = args.chet
        self.snp = args.snp

    def run(self):
        dirList = [self.Dir + '/' + x for x in os.listdir(self.Dir)]
        for subdir in dirList:
            if not os.path.isdir(subdir):
                continue
            subList = [subdir + '/' + x for x in os.listdir(subdir)]
            self.OnePed(subList)
        self.FlushBamout()
        self.FlushIGV()

    def OnePed(self, subList):
        PedFil = None
        PSAPFil = None
        for _file in subList:
            if _file.endswith('.ped'):
                PedFil = _file
            if _file.endswith('.csv'):
                PSAPFil = _file
        if PedFil == None or PSAPFil == None:
            return None
        print PedFil.split('/')[-1], PSAPFil.split('/')[-1]
        self.ReadCSV(PSAPFil, PedFil)

    def ReadCSV(self, PSAPFil, ped):
        pedigree = Pedigree(ped)
        Proband, Father, Mother = pedigree.GuessProband()
        ProbandBam = self.bampath.search(Proband)
        FatherBam = self.bampath.search(Father)
        MotherBam = self.bampath.search(Mother)
        res = []
        with open(PSAPFil, 'rb') as csvfile:
            reader = csv.reader(csvfile)
            header = reader.next()
            idx_Chr = header.index("Chrom")
            idx_Pos = header.index("Pos")
            idx_Ref = header.index("Ref")
            idx_Alt = header.index("Alt")
            idx_Proband = header.index(Proband)
            idx_DModel = header.index("Dz.Model." + Proband)
            idx_PopScore = header.index("popScore." + Proband)
            for row in reader:
                var = Variant(row, idx_Chr, idx_Pos, idx_Ref,
                              idx_Alt, idx_Proband, idx_DModel, idx_PopScore)
                if var.needIGV(self.popScore, self.snp, self.chet):
                    res.append(var.Out2L())
                    self.VarList2_hand.write("{}\t{},{},{},{},{},{}\n".format(var.Out2I(
                    ), ProbandBam.BamName, ProbandBam.BamOutName, FatherBam.BamName, FatherBam.BamOutName, MotherBam.BamName, MotherBam.BamOutName))
            self.VarList1_hand.write("{}\t{}\n".format(
                ProbandBam.FullPath, '\t'.join(res)))
            self.VarList1_hand.write("{}\t{}\n".format(
                FatherBam.FullPath, '\t'.join(res)))
            self.VarList1_hand.write("{}\t{}\n".format(
                MotherBam.FullPath, '\t'.join(res)))
        return res

    def FlushBamout(self):
        fout = self.Script1
        fout.write('REF={}\n'.format(EXM_REF))
        fout.write('CMD={}\n'.format(BAMOUT_CMD))
        fout.write('InpFil={}\n'.format(self.VarList1))
        fout.write("NumJob=`wc -l $InpFil|cut -f1 -d ' '`\n\n")
        fout.write("mkdir -p Bamouts;cd Bamouts\n")
        fout.write(
            "seq $NumJob| parallel -j {} --eta $CMD -i $InpFil -a {} -r $REF\n\n".format(self.parallel, "{}"))
        fout.write("cd ../;find `pwd`/Bamouts/ -name \"*.bamout.bam\" > Bamout.bam.list\n")
        fout.write("cat {} {} > ALL.bam.list\n".format(
            self.bampath.bamlocationfile, "Bamout.bam.list"))
        fout.close()

    def FlushIGV(self):
        fout = self.Script2
        fout.write(
            '{} -b ALL.bam.list -v {}'.format(IGV_GENERATE_CMD, self.VarList2))
        fout.close()


class Indv:
    def __init__(self, FamID, SampleID, FatherID, MotherID, Sex, Affected):
        self.FamID = FamID
        self.SampleID = SampleID
        self.FatherID = FatherID
        self.MotherID = MotherID
        self.Sex = Sex
        self.Affected = Affected

    def show(self):
        print '\t'.join([self.FamID, self.SampleID, self.FatherID, self.MotherID, self.Sex, self.Affected])

# Type: Singleton, Trios and complex Fam.


class Pedigree:
    def __init__(self, pedfile):
        fin = open(pedfile, 'rb')
        self.Indvs = []
        for l in fin:
            if l.startswith('#'):
                continue
            llist = l.strip().split('\t')
            FamID, SampleID, FatherID, MotherID, Sex, Affected = llist[:6]
            indv = Indv(FamID, SampleID, FatherID, MotherID, Sex, Affected)
            self.Indvs.append(indv)
        self.FamID = FamID

    def GuessProband(self):  # Proband should be affected and have max number of parents
        ProbandCandidates = []
        for Indv in self.Indvs:
            # Indv.show()
            if Indv.Affected == '2':
                NumParents = self.GetNumParents(Indv)
                ProbandCandidates.append((Indv, NumParents))
                # return Indv.SampleID
        if len(ProbandCandidates) == 0:
            print "Warning: Can't find a Proband at FAM: {}".format(self.FamID)
        else:
            Proband = sorted(ProbandCandidates,
                             key=lambda x: x[1], reverse=True)[0][0]
            return Proband.SampleID, Proband.FatherID, Proband.MotherID

        return self.Indvs[0].SampleID

    def GetNumParents(self, Sample):
        F = 1 if str(Sample.FatherID) != "0" else 0
        M = 1 if str(Sample.MotherID) != "0" else 0
        return sum([F, M])


class BAM:
    def __init__(self, fpath):
        self.FullPath = fpath.strip()
        self.BamName = self.FullPath.split('/')[-1]
        self.BamOutName = '.'.join(
            self.BamName.split('.')[:-1]) + ".bamout.bam"
        try:
            self.Sample = SampleName.search(self.BamName).group(1)
        except:
            print "SampleName {} Not get matched.".format(self.BamName)
            exit()
        #self.Sample = self.BamName


class BamLocation:
    def __init__(self, bamlocationfile):
        fin = open(bamlocationfile, 'rb')
        self.bamlocationfile = os.path.abspath(bamlocationfile)
        self.Bams = {}
        for l in fin:
            if 'bamout' in l:
                continue
            bam = BAM(l)
            self.Bams[bam.Sample] = bam

    def search(self, sample):
        return self.Bams.get(sample, 'NA')


class Variant:
    def __init__(self, row, idx_Chr, idx_Pos, idx_Ref, idx_Alt, idx_Proband, idx_DModel, idx_PopScore):
        self.Chrom = row[idx_Chr]
        self.Pos = int(row[idx_Pos])
        self.Ref = row[idx_Ref]
        self.Alt = row[idx_Alt].split(',')
        self.Sample = row[idx_Proband]
        self.DModel = row[idx_DModel]
        self.PopScore = row[idx_PopScore]
        self.start = max(0, self.Pos - 150)
        self.end = self.Pos + 150
    # INDELs and Chet with top PopScore needs IGV comfirm.

    def needIGV(self, popscoreCutoff, snp, chet):
        if float(self.PopScore) <= popscoreCutoff:
            if snp or (len(self.Alt)) > 1 or (len(self.Ref) != len(self.Alt[0])):
                return True
            elif chet and (self.DModel == 'REC-chet'):
                return True
        else:
            return False

    def Out2L(self):  # To run_Bamout.sh
        return "-L {}:{}-{}".format(self.Chrom, self.start, self.end)

    def Out2I(self):  # To run_IGV.sh
        return "{}\t{}".format(self.Chrom, self.Pos)


def main():
    args = GetOptions()
    instance = PSAP_PrepareIGVInput(args)
    instance.run()
    return


if __name__ == '__main__':
    main()
