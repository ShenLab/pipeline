#!/usr/bin/env python

# The purpose of this script is to count genotypes at each variant
#    -o/--out    <required>    The use should specify a base name for the output files.The script outputs the filtered results as both a tab delimited file and a vcf file. The script also produces a log file and an auxilary file that contains the results of every test of every variant as "TRUE" or "FALSE" (primarily for debugging purposes).

import sys
import os
from optparse import OptionParser
parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",
                  help="Input VCF file", metavar="VCFfile")
parser.add_option("-o", "--out", dest="OutputFileName",
                  help="Base of names for output files", metavar="OutputFileName")

(options, args) = parser.parse_args()


# Assign input and output files
VCF = open(options.VCFfile, 'r')
BaseName = str(options.OutputFileName)
TabOutputFilename = BaseName + '.counts.tsv'
Output = open(TabOutputFilename, 'w')

# Define sample names using user-defined parameters

# start log file

# Start output table
headerlist = ['Chromosome', 'Position', 'ID', 'REF', 'ALT', 'Gene', 'VariantFunction', 'VariantClass', 'AAchange', '1KG.score', 'ESP.score', 'SIFTscore', 'SIFTprediction',
              'PP2score', 'PP2prediction', 'MTscore', 'MTprediction', 'GERPscore', 'PHYLOPscore', 'CADDscore', 'QD', 'QUAL', 'INFO', 'Ref', 'Het', 'Alt', 'NoCall', 'SampleNum']
Output.write("\t".join(headerlist) + "\n")

ChrCount = 0
VarCount = 0
# Read VCF file
for line in VCF:
    # Map column name to number, and then find column numbers of each set of
    # trios
    if '#' not in line:
        # Variant must first pass 1KG and GO-ESP frequencies, MQ0 threshold,
        # and be exonic
        linelist = line.split("\t")
        ChrPresent = linelist[0]
        if ChrCount != ChrPresent:
            print
            sys.stdout.write('Chromosome ' + str(ChrPresent) + ' ')
            ChrCount = ChrPresent
            VarCount = 0
        VarCount = VarCount + 1
        if VarCount == 1:
            sys.stdout.write('.')
            sys.stdout.flush()
        if VarCount == 1000:
            VarCount = 0
        QUAL = linelist[5]
        INFOstring = linelist[7]
        INFOcolumnList = INFOstring.split(";")
        INFOdict = {}
        for element in INFOcolumnList:
            if '=' in element:
                FieldName, FieldValue = element.split('=', 1)
                INFOdict[FieldName] = FieldValue

        # Get variant data
        QDnumber = float(INFOdict.get('QD', 0))
        DPnumber = float(INFOdict.get('DP', 0))
        MQ0number = float(INFOdict.get('MQ0', 0))
        GeneName = INFOdict.get('GeneName', 'NA')
        VariantFunction = INFOdict.get('VarFunc', 'none')
        VariantFunctionList = VariantFunction.split(',')
        VariantClass = INFOdict.get('VarClass', 'none')
        VariantClassList = VariantClass.split(',')
        AAchange = INFOdict.get('AAChange', 'NA')
        KGscore = float(INFOdict.get('1KGfreq', 0))
        ESPscore = float(INFOdict.get('ESPfreq', 0))
        SIFTscore = float(INFOdict.get('SIFTscr', 0))
        SIFTprediction = INFOdict.get('SIFTprd', 'NA')
        PP2score = float(INFOdict.get('PP2scr', 0))
        PP2prediction = INFOdict.get('PP2prd', 'NA')
        MTscore = float(INFOdict.get('MutTscr', 0))
        MTprediction = INFOdict.get('MutTprd', 'NA')
        GERPscore = float(INFOdict.get('GERP', 0))
        PHYLOPscore = float(INFOdict.get('PhyloP', 0))
        CADDscore = float(INFOdict.get('CADD', 0))

        # Get number of alternate alleles
        AltAlls = linelist[4]
        AltAlls = AltAlls.split(",")
        AltNum = len(AltAlls)

        # Get Getnotypes(
        QualityString = [linelist[i].strip() for i in range(9, len(linelist))]
        QualityList = [i.split(':') for i in QualityString]
        GTList = [QualityList[i][0] for i in range(0, len(QualityString))]

        # Count GTs
        SampleNum = len(linelist) - 9
        RefCnt = 0
        NclCnt = 0
        RefCnt = GTList.count('0/0')
        NclCnt = GTList.count('./.')

        # Get number of alternate alleles
        AltAlls = linelist[4]
        AltAlls = AltAlls.split(",")
        AltNum = len(AltAlls)

        if AltNum > 1:
            for refnum in range(0, AltNum):
                HetCnt = 0
                AltCnt = 0
                HetStr = '0/' + str(refnum + 1)
                AltStr = str(refnum + 1) + '/' + str(refnum + 1)
                HetCnt = GTList.count(HetStr)
                AltCnt = GTList.count(AltStr)

                SttList = linelist[0:5]
                SttList[4] = AltAlls[refnum]
                OutputList = SttList + [GeneName, VariantFunction, VariantClass, AAchange, KGscore, ESPscore, SIFTscore, SIFTprediction, PP2score, PP2prediction,
                                        MTscore, MTprediction, GERPscore, PHYLOPscore, CADDscore, QDnumber, QUAL] + linelist[7:8] + [RefCnt, HetCnt, AltCnt, NclCnt, SampleNum]
                OutputList = [str(i) for i in OutputList]
                OutputString = "\t".join(OutputList)
                Output.write(OutputString + "\n")

        else:
            HetCnt = 0
            AltCnt = 0

            HetCnt = GTList.count('0/1')
            AltCnt = GTList.count('1/1')

            OutputList = linelist[0:5] + [GeneName, VariantFunction, VariantClass, AAchange, KGscore, ESPscore, SIFTscore, SIFTprediction, PP2score, PP2prediction,
                                          MTscore, MTprediction, GERPscore, PHYLOPscore, CADDscore, QDnumber, QUAL] + linelist[7:8] + [RefCnt, HetCnt, AltCnt, NclCnt, SampleNum]
            OutputList = [str(i) for i in OutputList]
            OutputString = "\t".join(OutputList)
            Output.write(OutputString + "\n")
