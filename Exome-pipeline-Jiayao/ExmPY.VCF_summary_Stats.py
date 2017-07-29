#!/usr/bin/env python


# this script takes as its input a vcf file and outputs a tab delimited
# containing various stats for each sample in the vcf, e.g. number of
# SNVs, number of InDels Ti/Tv ratios etc.
import gzip
from optparse import OptionParser
parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",
                  help="input VCF file", metavar="VCFfile")
parser.add_option("-o", "--output", dest="OutputFileName",
                  help="user specified name of file for output to be written", metavar="OutputFileName")


parser.add_option("-F", "--filter", action='store_true',
                  dest="FilterVars", help="Only count good quality variants")

(options, args) = parser.parse_args()

OutFile = str(options.OutputFileName)
if not OutFile.endswith('.tsv'):
    OutFile += '.stats.tsv'
Output = open(OutFile, 'w')
# MigFil=open('TESTLog.tsv','w')
#MigHeader=['Sample','CHR', 'Pos', 'REF', 'ALT', 'GT', 'Exon', 'Class', 'Heterozygous', 'Homozygous' 'UnCalled', 'Indel', 'SNV', 'INFO']
# MigFil.write("\t".join(MigHeader)+"\n")


FilterVariants = options.FilterVars

FilTyp = str(options.VCFfile)
FilTyp = FilTyp.split(".")
FilTyp = str(FilTyp[-1])
if FilTyp != 'gz'and FilTyp != 'vcf':
    print "Incorrect input file type"
    sys.exit(1)

Nucleotides = ['A', 'C', 'G', 'T']
MutUnknown = ['unknown']
MutSilent = ['synonymousSNV']
MutMissense = ['nonsynonymousSNV']
MutNonsense = ['stopgain', 'stoploss']
MutNonframeshift = ['nonframeshiftinsertion', 'nonframeshiftdeletion']
MutFrameshift = ['frameshiftinsertion', 'frameshiftdeletion']
CodingCodes = ['splicing', 'exonic', 'exonic;splicing']
SplicingCodes = ['splicing', 'exonic;splicing']
CodingTypes = ['Coding variants', 'Non-coding Variants', 'All variants']
BadFilters = ['QD_Bad_SNP', 'FS_Bad_SNP',
              'FS_Mid_SNP;QD_Mid_SNP', 'LowQD_Indel']

print "Running VCF stats..."

WildType = False
Heterozygous = False
Homozygous = False
UnCalled = False
Indel = False
SNV = False
Transversion = False
Transition = False
Silent = False
Missense = False
Nonsense = False
Unknown = False
Known = False
Novel = False
FrameShift = False
Rare = False

for WhichCode in CodingTypes:
    print WhichCode
    if FilTyp == 'gz':
        VCF = gzip.open(options.VCFfile, 'r')
    elif FilTyp == 'vcf':
        VCF = open(options.VCFfile, 'r')
    # Read VCF file
    for line in VCF:
        line = line.strip()
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
            KnownTiCount = [0] * SamLength
            KnownTvCount = [0] * SamLength
            NovelTiCount = [0] * SamLength
            NovelTvCount = [0] * SamLength
            SilentCount = [0] * SamLength
            MissenseCount = [0] * SamLength
            NonsenseCount = [0] * SamLength
            UnknownCount = [0] * SamLength
            FrameShiftCount = [0] * SamLength
            WildTypeCount = [0] * SamLength
            HomozygousCount = [0] * SamLength
            HeterozygousCount = [0] * SamLength
            RareCount = [0] * SamLength
            RareSilentCount = [0] * SamLength
            KnownTiTvRat = [0] * SamLength
            NovelTiTvRat = [0] * SamLength
            TotalTiTvRat = [0] * SamLength
            UnCalledCount = [0] * SamLength
        # Start Countin
        if '#' not in line:
            # Variant must first pass 1KG and GO-ESP frequencies, MQ0
            # threshold, and be exonic
            linelist = line.split("\t")
            VariantFilter = linelist[6]
            INFOstring = linelist[7]
            INFOcolumnList = INFOstring.split(";")
            INFOdict = {}
            for element in INFOcolumnList:
                if '=' in element:
                    FieldName, FieldValue = element.split('=', 1)
                    INFOdict[FieldName] = FieldValue

            # Get values for later
            QDnumber = float(INFOdict.get('QD', 0))

            ExAClist = str(INFOdict.get('ExACfreq', 2))
            ExAClist = ExAClist.split(",")

            MutationFunct = str(INFOdict.get('VarFunc', 'none'))
            MutationFunct = MutationFunct.split(",")
            MutationFunct = MutationFunct[0]

            MutationClass = str(INFOdict.get('VarClass', 'none'))
            MutationClass = MutationClass.split(",")

            #Known or Novel
            ID = str(linelist[2])

            # check whether it is a coding variant and then run only if
            # appropriate
            if \
                    WhichCode == "Coding variants" and MutationFunct in CodingCodes or \
                    WhichCode == "Non-coding Variants" and MutationFunct not in CodingCodes or \
                    WhichCode == "All variants":
                # first run a dummy sample to count the "total" number of
                # variants, this is actually just the number of loci and
                # ignores mulit-allelic factors, then loop through each sample
                # column
                for SampleNumber in range(0, SamLength):

                    # first get the genotype for the sample
                    SampleColumn = SampleNumber + 8
                    if SampleNumber == 0:
                        # this is for the total so make the genotype the first
                        # alternate allele
                        GT = "0/1"
                    else:
                        InfoList = linelist[SampleColumn].strip()
                        InfoListSPL = InfoList.split(":")
                        GT = InfoListSPL[0]
                        if str(GT) == "./.":
                            # the genotype was not called...
                            UnCalled = True
                        if str(GT) == ".":
                            # the genotype was not called...
                            UnCalled = True
                        if str(GT) == "0/0":
                            # the genotype was not called...
                            WildType = True

                    # continue only if this is a variant (not 0/0 or uncalled)
                    # otherwise move on to the next sample
                    if not WildType and not UnCalled:
                        if 'DB' in INFOcolumnList or ID != ".":
                            Known = True
                        else:
                            Novel = True

                        AllList = set(GT.split("/"))
                        if len(AllList) == 2:
                            Heterozygous = True
                        else:
                            Homozygous = True
                        # Get the alternate allele. Only use the first
                        # alternate allele in cases where the sample has two
                        # different variants (e.g. 1/2)
                        AltAllele = GT.split("/")
                        if "0" in AltAllele:
                            AltAllele.remove("0")
                        AltAllele = AltAllele[0]
                        AlNum = int(AltAllele) - 1

                        #Indel or SNP
                        MutNum = min(AlNum, len(MutationClass) - 1)
                        if MutationClass[MutNum] in MutFrameshift or MutationClass[MutNum] in MutNonframeshift:
                            Indel = True
                        else:
                            SNV = True

                        #Ti or Tv
                        REF = str(linelist[3])
                        ALT = linelist[4].split(",")
                        ALT = ALT[AlNum]
                        if (REF[0] == 'A' and ALT[0] == 'G') or (REF[0] == 'G' and ALT[0] == 'A') or (REF[0] == 'C' and ALT[0] == 'T') or (REF[0] == 'T' and ALT[0] == 'C'):
                            Transition = True
                        if (REF[0] == 'C' and ALT[0] == 'A') or (REF[0] == 'A' and ALT[0] == 'C') or (REF[0] == 'G' and ALT[0] == 'C') or (REF[0] == 'C' and ALT[0] == 'G') or (REF[0] == 'T' and ALT[0] == 'A') or (REF[0] == 'A' and ALT[0] == 'T') or (REF[0] == 'T' and ALT[0] == 'G') or (REF[0] == 'G' and ALT[0] == 'T'):
                            Transversion = True

                        # Class
                        if MutationClass[MutNum] in MutSilent:
                            Silent = True
                        if MutationClass[MutNum] in MutMissense:
                            Missense = True
                        if MutationClass[MutNum] in MutNonsense or MutationFunct in SplicingCodes:
                            Nonsense = True
                        if MutationClass[MutNum] in MutUnknown:
                            Unknown = True
                        if MutationClass[MutNum] in MutFrameshift:
                            FrameShift = True

                        # Frequency
                        ExACscore = ExAClist[min(AlNum, len(ExAClist) - 1)]
                        if str(ExACscore) == ".":
                            ExACscore = 0.005
                        ExACscore = float(ExACscore)
                        if ExACscore <= 0.01:
                            Rare = True
                        # MigList=[Samples[SampleNumber-1]]+linelist[0:2]+[REF,ALT,GT,MutationFunct,MutationClass[MutNum],str(Heterozygous),str(Homozygous),str(UnCalled),str(Indel),str(SNV),INFOstring]
                        # MigFil.write("\t".join(#MigList)+"\n")

                    # add to counts if filter pass...
                    if all(str(i) not in VariantFilter for i in BadFilters) or not FilterVariants:
                        if SNV:
                            SNVcount[SampleNumber] = SNVcount[SampleNumber] + 1
                        if Indel:
                            InDelCount[SampleNumber] = InDelCount[SampleNumber] + 1
                        if Known:
                            KnownCount[SampleNumber] = KnownCount[SampleNumber] + 1
                        if Novel:
                            NovelCount[SampleNumber] = NovelCount[SampleNumber] + 1
                        if Known and Transition:
                            KnownTiCount[SampleNumber] = KnownTiCount[SampleNumber] + 1
                        if Known and Transversion:
                            KnownTvCount[SampleNumber] = KnownTvCount[SampleNumber] + 1
                        if Novel and Transition:
                            NovelTiCount[SampleNumber] = NovelTiCount[SampleNumber] + 1
                        if Novel and Transversion:
                            NovelTvCount[SampleNumber] = NovelTvCount[SampleNumber] + 1
                        if Silent:
                            SilentCount[SampleNumber] = SilentCount[SampleNumber] + 1
                        if Missense:
                            MissenseCount[SampleNumber] = MissenseCount[SampleNumber] + 1
                        if Nonsense:
                            NonsenseCount[SampleNumber] = NonsenseCount[SampleNumber] + 1
                        if Unknown:
                            UnknownCount[SampleNumber] = UnknownCount[SampleNumber] + 1
                        if FrameShift:
                            FrameShiftCount[SampleNumber] = FrameShiftCount[SampleNumber] + 1
                        if Homozygous:
                            HomozygousCount[SampleNumber] = HomozygousCount[SampleNumber] + 1
                        if Heterozygous:
                            HeterozygousCount[SampleNumber] = HeterozygousCount[SampleNumber] + 1
                        if WildType:
                            WildTypeCount[SampleNumber] = WildTypeCount[SampleNumber] + 1
                        if UnCalled:
                            UnCalledCount[SampleNumber] = UnCalledCount[SampleNumber] + 1
                        if Rare:
                            RareCount[SampleNumber] = RareCount[SampleNumber] + 1
                        if Silent and Rare:
                            RareSilentCount[SampleNumber] = RareSilentCount[SampleNumber] + 1
                        SNV = False
                        Indel = False
                        Transversion = False
                        Transition = False
                        Known = False
                        Novel = False
                        Silent = False
                        Missense = False
                        Nonsense = False
                        FrameShift = False
                        Unknown = False
                        Homozygous = False
                        Heterozygous = False
                        WildType = False
                        UnCalled = False
                        Rare = False

    for i in range(0, SamLength):
        if int(KnownTvCount[i]) > 0:
            KnownTiTvRat[i] = float(KnownTiCount[i]) / float(KnownTvCount[i])
        if int(NovelTvCount[i]) > 0:
            NovelTiTvRat[i] = float(NovelTiCount[i]) / float(NovelTvCount[i])
        TotalTvCount = float(NovelTvCount[i]) + float(KnownTvCount[i])
        TotalTiCount = float(NovelTiCount[i]) + float(KnownTiCount[i])
        if TotalTvCount > 0:
            TotalTiTvRat[i] = float(TotalTiCount) / float(TotalTvCount)
    # MigFil.write("\n")
    VCF.close()
    # write output
    Output.write("\t".join([str(WhichCode)]) + "\n")
    Output.write("\t".join(['Samples:'] + [str(i) for i in Samples]) + "\n")
    Output.write("\t".join(['SNVs'] + [str(i) for i in SNVcount]) + "\n")
    Output.write("\t".join(['InDels'] + [str(i) for i in InDelCount]) + "\n")
    Output.write("\t".join(['Known'] + [str(i) for i in KnownCount]) + "\n")
    Output.write("\t".join(['Novel'] + [str(i) for i in NovelCount]) + "\n")
    Output.write("\t".join(['Known Ti'] + [str(i)
                                           for i in KnownTiCount]) + "\n")
    Output.write("\t".join(['Known TV'] + [str(i)
                                           for i in KnownTvCount]) + "\n")
    Output.write("\t".join(['Known Ti/TV ratio'] + [str(i)
                                                    for i in KnownTiTvRat]) + "\n")
    Output.write("\t".join(['Novel Ti'] + [str(i)
                                           for i in NovelTiCount]) + "\n")
    Output.write("\t".join(['Novel TV'] + [str(i)
                                           for i in NovelTvCount]) + "\n")
    Output.write("\t".join(['Novel Ti/TV ratio'] + [str(i)
                                                    for i in NovelTiTvRat]) + "\n")
    Output.write("\t".join(['Overall Ti/TV ratio'] + [str(i)
                                                      for i in TotalTiTvRat]) + "\n")
    if str(CodingTypes) != 'Non-coding Variants':
        Output.write("\t".join(['Silent'] + [str(i)
                                             for i in SilentCount]) + "\n")
        Output.write("\t".join(['Missense'] + [str(i)
                                               for i in MissenseCount]) + "\n")
        Output.write("\t".join(['Nonsense'] + [str(i)
                                               for i in NonsenseCount]) + "\n")
        Output.write("\t".join(['UnKnown'] + [str(i)
                                              for i in UnknownCount]) + "\n")
        Output.write("\t".join(['Frameshift'] + [str(i)
                                                 for i in FrameShiftCount]) + "\n")
    Output.write("\t".join(['Homozygous'] + [str(i)
                                             for i in HomozygousCount]) + "\n")
    Output.write("\t".join(['Heterozygous'] + [str(i)
                                               for i in HeterozygousCount]) + "\n")
    Output.write("\t".join(['Reference'] + [str(i)
                                            for i in WildTypeCount]) + "\n")
    Output.write("\t".join(['Not Called'] + [str(i)
                                             for i in UnCalledCount]) + "\n")
    Output.write("\t".join(
        ['Rare (AAF<0.01 or missing in ExAC)'] + [str(i) for i in RareCount]) + "\n")
    if str(CodingTypes) != 'Non-coding Variants':
        Output.write("\t".join(['Rare Silent'] + [str(i)
                                                  for i in RareSilentCount]) + "\n")
    Output.write("\t\n")
