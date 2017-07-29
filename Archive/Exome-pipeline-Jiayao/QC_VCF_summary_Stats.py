#!/usr/bin/env python
from __future__ import division

import gzip
import matplotlib as mpl
mpl.use('Agg')
import os
from optparse import OptionParser
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import numpy as np

'''
script to summarize vcf files

Input
----------
:vcf

Output
----------
summary: variants counts summary 
pdf: visulized summary 

'''


parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",
                  help="input VCF file", metavar="VCFfile")
parser.add_option("-q", "--qc_folder", dest="QC_folder",
                  help="user specified name of folder for output to be written", metavar="QC_folder")

(options, args) = parser.parse_args()
QC_folder = os.path.abspath(options.QC_folder)

outname = os.path.join(QC_folder, 'sample.summary')

Output = open(outname, 'w')

Nucleotides = ['A', 'C', 'G', 'T']
MutUnknown = ['unknown']
MutSilent = ['synonymous_SNV']
MutMissense = ['nonsynonymous_SNV']
MutNonsense = ['stopgain', 'stoploss']
MutNonframeshift = ['nonframeshift_insertion', 'nonframeshift_deletion']
MutFrameshift = ['frameshift_insertion', 'frameshift_deletion']
CodingCodes = ['splicing', 'exonic', 'exonic-splicing']
CodingTypes = ['Coding variants', 'Non-coding Variants', 'All variants']

print '-' * 50
print 'QC summary started'

for WhichCode in CodingTypes:
    print WhichCode
    if options.VCFfile.endswith('vcf.gz'):
        VCF = gzip.open(options.VCFfile, 'r')
    else:
        VCF = open(options.VCFfile, 'r')
    # Read VCF file
    for line in VCF:
        # Map column name to number, and then find column numbers of each set
        # of trios
        if '#CHROM' in line:
            Samples = line.split("\t")
            Samples = Samples[8:]
            Samples[0] = 'Total'
            SamLength = len(Samples)
            SNVcount = [0] * SamLength
            InDelCount = [0] * SamLength
            KnownCount = [0] * SamLength
            NovelCount = [0] * SamLength
            knownTiCount = [0] * SamLength
            knownTvCount = [0] * SamLength
            novelTiCount = [0] * SamLength
            novelTvCount = [0] * SamLength
            silentCount = [0] * SamLength
            missenseCount = [0] * SamLength
            nonsenseCount = [0] * SamLength
            unknownCount = [0] * SamLength
            frameShiftCount = [0] * SamLength
            homozygousCount = [0] * SamLength
            heterozygousCount = [0] * SamLength
            referenceCount = [0] * SamLength
            KGRareCount = [0] * SamLength
            ESPRareCount = [0] * SamLength
            RareCount = [0] * SamLength
            ExACRareCount = [0] * SamLength
            NoAAF = [0] * SamLength
            KnownTiTvRat = [0] * SamLength
            NovelTiTvRat = [0] * SamLength
            TotalTiTvRat = [0] * SamLength
            UnCalled = [0] * SamLength
        # Start Countin
        if '#' not in line:
            # Variant must first pass 1KG and GO-ESP frequencies, MQ0
            # threshold, and be exonic
            linelist = line.split("\t")
            INFOstring = linelist[7]
            INFOcolumnList = INFOstring.split(";")
            INFOdict = {}
            for element in INFOcolumnList:
                if '=' in element:
                    FieldName, FieldValue = element.split('=', 1)
                    if FieldName not in INFOdict:
                        INFOdict[FieldName] = FieldValue
                    else:
                        INFOdict[FieldName] = INFOdict[FieldName] + \
                            ',' + FieldValue

            # Get values for later
            KGseq = map(float, [rate for rate in INFOdict.get(
                '1000g2015aug_all', '2').split(',') if rate != '.'])
            ESPseq = map(float, [rate for rate in INFOdict.get(
                'esp6500siv2_all', '2').split(',') if rate != '.'])
            ExACseq = map(float, [rate for rate in INFOdict.get(
                'ExAC_ALL', '2').split(',') if rate != '.'])
            KGseq.append(0)
            ESPseq.append(0)
            ExACseq.append(0)

            KGscore = max(KGseq)
            ESPscore = max(ESPseq)
            ExACscore = max(ExACseq)

            MutationFunct = str(INFOdict.get(
                'Func.refGene', 'none').split(',')[0])
            MutationClass = str(INFOdict.get(
                'ExonicFunc.refGene', 'none').split(',')[0])
            ID = str(linelist[2])
            REF = linelist[3].split(",")
            REF = [str(i) for i in REF]
            ALT = linelist[4].split(",")
            ALT = [str(i) for i in ALT]
            #Known or novel
            Known = False
            if INFOdict['avsnp147'].split(',')[0] != '.':
                Known = True
            if ID != ".":
                Known = True
            #Snp or Indel
            InDel = True
            if all(i in Nucleotides for i in REF) and all(i in Nucleotides for i in ALT):
                InDel = False
            #Ti or Tv
            Transversion = False
            Transition = False
            if (REF[0] == 'A' and ALT[0] == 'G') or (REF[0] == 'G' and ALT[0] == 'A') or (REF[0] == 'C' and ALT[0] == 'T') or (REF[0] == 'T' and ALT[0] == 'C'):
                Transition = True
            if (REF[0] == 'C' and ALT[0] == 'A') or (REF[0] == 'A' and ALT[0] == 'C') or (REF[0] == 'G' and ALT[0] == 'C') or (REF[0] == 'C' and ALT[0] == 'G') or (REF[0] == 'T' and ALT[0] == 'A') or (REF[0] == 'A' and ALT[0] == 'T') or (REF[0] == 'T' and ALT[0] == 'G') or (REF[0] == 'G' and ALT[0] == 'T'):
                Transversion = True
            codingPass = False
            if WhichCode == "Coding variants" and MutationFunct in CodingCodes:
                codingPass = True
            if WhichCode == "Non-coding Variants" and MutationFunct not in CodingCodes:
                codingPass = True
            if WhichCode == "All variants":
                codingPass = True

            if codingPass:
                for i in range(0, SamLength):
                    ColNum = i + 8
                    if i != 0:
                        InfoList = linelist[ColNum].strip()
                        if "./." not in InfoList:
                            InfoListSPL = InfoList.split(":")
                            GT = InfoListSPL[0]
                            if str(GT) != "0/0":
                                AllList = GT.split("/")
                                if "0" in AllList:
                                    heterozygousCount[i] = heterozygousCount[i] + 1
                                else:
                                    homozygousCount[i] = homozygousCount[i] + 1
                        else:
                            UnCalled[i] = UnCalled[i] + 1
                    else:
                        GT = "0/0"
                    if str(GT) != "0/0" or i == 0:
                        if InDel:
                            InDelCount[i] = InDelCount[i] + 1
                        else:
                            SNVcount[i] = SNVcount[i] + 1
                        if Known:
                            KnownCount[i] = KnownCount[i] + 1
                            if Transition:
                                knownTiCount[i] = knownTiCount[i] + 1
                            if Transversion:
                                knownTvCount[i] = knownTvCount[i] + 1
                        else:
                            NovelCount[i] = NovelCount[i] + 1
                            if Transition:
                                novelTiCount[i] = novelTiCount[i] + 1
                            if Transversion:
                                novelTvCount[i] = novelTvCount[i] + 1
                        if MutationClass in MutSilent:
                            silentCount[i] = silentCount[i] + 1
                        if MutationClass in MutMissense:
                            missenseCount[i] = missenseCount[i] + 1
                        if MutationClass in MutNonsense:
                            nonsenseCount[i] = nonsenseCount[i] + 1
                        if MutationClass in MutUnknown:
                            unknownCount[i] = unknownCount[i] + 1
                        if MutationClass in MutFrameshift:
                            frameShiftCount[i] = frameShiftCount[i] + 1
                        if KGscore <= 0.01:
                            KGRareCount[i] = KGRareCount[i] + 1
                        if ESPscore <= 0.01:
                            ESPRareCount[i] = ESPRareCount[i] + 1
                        if ExACscore <= 0.01:
                            ExACRareCount[i] = ExACRareCount[i] + 1

                        if KGscore <= 0.01 or ESPscore <= 0.01:
                            RareCount[i] = RareCount[i] + 1
                        if KGscore == 2 and ESPscore == 2:
                            NoAAF[i] = NoAAF[i] + 1
    VCF.close()
    for i in range(0, SamLength):
        if int(knownTvCount[i]) > 0:
            KnownTiTvRat[i] = float(knownTiCount[i]) / float(knownTvCount[i])
        if int(novelTvCount[i]) > 0:
            NovelTiTvRat[i] = float(novelTiCount[i]) / float(novelTvCount[i])
        TotalTvCount = float(novelTvCount[i]) + float(knownTvCount[i])
        TotalTiCount = float(novelTiCount[i]) + float(knownTiCount[i])
        if TotalTvCount > 0:
            TotalTiTvRat[i] = float(TotalTiCount) / float(TotalTvCount)
    # write output
    Output.write("\t".join([str(WhichCode)]) + "\n")
    Output.write("\t".join(['Samples:'] + [str(i) for i in Samples]))
    Output.write("\t".join(['SNVs'] + [str(i) for i in SNVcount]) + "\n")
    Output.write("\t".join(['InDels'] + [str(i) for i in InDelCount]) + "\n")
    Output.write("\t".join(['Known'] + [str(i) for i in KnownCount]) + "\n")
    Output.write("\t".join(['Novel'] + [str(i) for i in NovelCount]) + "\n")
    Output.write("\t".join(['known Ti'] + [str(i)
                                           for i in knownTiCount]) + "\n")
    Output.write("\t".join(['known TV'] + [str(i)
                                           for i in knownTvCount]) + "\n")
    Output.write("\t".join(['known Ti/TV ratio'] + [str(i)
                                                    for i in KnownTiTvRat]) + "\n")
    Output.write("\t".join(['novel Ti'] + [str(i)
                                           for i in novelTiCount]) + "\n")
    Output.write("\t".join(['novel TV'] + [str(i)
                                           for i in novelTvCount]) + "\n")
    Output.write("\t".join(['novel Ti/TV ratio'] + [str(i)
                                                    for i in NovelTiTvRat]) + "\n")
    Output.write("\t".join(['total Ti/TV ratio'] + [str(i)
                                                    for i in TotalTiTvRat]) + "\n")
    if str(CodingTypes) != 'Non-coding Variants':
        Output.write("\t".join(['Silent'] + [str(i)
                                             for i in silentCount]) + "\n")
        Output.write("\t".join(['Missense'] + [str(i)
                                               for i in missenseCount]) + "\n")
        Output.write("\t".join(['Nonsense'] + [str(i)
                                               for i in nonsenseCount]) + "\n")
        Output.write("\t".join(['Unknown'] + [str(i)
                                              for i in unknownCount]) + "\n")
        Output.write("\t".join(['Frameshift'] + [str(i)
                                                 for i in frameShiftCount]) + "\n")
    Output.write("\t".join(['Homozygous'] + [str(i)
                                             for i in homozygousCount]) + "\n")
    Output.write("\t".join(['Heterozygous'] + [str(i)
                                               for i in heterozygousCount]) + "\n")
    Output.write("\t".join(['Rare (AAF<0.01) - eKG'] +
                           [str(i) for i in KGRareCount]) + "\n")
    Output.write("\t".join(['Rare (AAF<0.01) - ExAC'] +
                           [str(i) for i in KGRareCount]) + "\n")
    Output.write("\t".join(['Rare (AAF<0.01) - ESP'] +
                           [str(i) for i in ESPRareCount]) + "\n")
    Output.write("\t".join(
        ['Rare (AAF<0.01) - 1KG or ESP'] + [str(i) for i in RareCount]) + "\n")
    Output.write("\t".join(['No AAF'] + [str(i) for i in NoAAF]) + "\n")
    Output.write("\t".join(['not called'] + [str(i) for i in UnCalled]) + "\n")
    Output.write("\t\n")
Output.close()


plot_caterogy = ['SNVs', 'InDels', 'known Ti/TV ratio',
                 'novel Ti/TV ratio', 'Rare (AAF<0.01) - ExAC']

print outname
with open(outname) as f:
    for line in f:
        line = line.strip()

        if line in CodingTypes:
            output_type = line
            continue
        if line == '':
            continue
        if line.startswith('Samples'):
            samples = line.split()
            print '# of samples in vcf:', len(samples) - 2
            continue

        info = line.split('\t')
        caterogy = info[0]
        if caterogy in plot_caterogy:
            print 'plot: {} - {}'.format(output_type, caterogy)

            hist = map(float, info[2:])
            pdfname = os.path.join(
                QC_folder, output_type + '-' + caterogy.replace('/', '') + '.pdf')
            pdf = PdfPages(pdfname)
            fig, ax = plt.subplots(dpi=100)
            rects1 = ax.hist(hist, bins=50)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.set_title(output_type + ' ' + caterogy + '\n' + 'mean:' + str(np.round(np.mean(hist), decimals=3)) + ' ' +
                         'median:' + str(np.round(np.median(hist), decimals=3)))
            plt.show()
            plt.tight_layout()

            pdf.savefig(bbox_inches='tight')
            pdf.close()
            plt.close()


print 'QC summary finished'
print '-' * 50
print ''
