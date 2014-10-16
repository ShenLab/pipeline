#!/usr/bin/env python
#$ -N CustFilt -cwd

# The purpose of this script is to filter a VCF file for variants based on user-specified genotypes for each sample.
# For example: To filter according to an autosomal dominant inheritance model, filter all mutations that are homozygous reference (0/0) in one parent and heterozygous (0/1) in the other parent and the proband.
#    -v/--vcf    <required>    The script requires an input vcf
#    -o/--out    <required>    The use should specify a base name for the output files.The script outputs the filtered results as both a tab delimited file and a vcf file. The script also produces a log file and an auxilary file that contains the results of every test of every variant as "TRUE" or "FALSE" (primarily for debugging purposes).
#    -a/--alt ...
#    -e/--het ...
#    -r/--ref ... <required>    For each genotype - homozygous alternate, heterozygous, homozygous reference, the user must supply the sample ids as in the vcf for each genotype. If there is more than one sample per a genotype they should be separated by a comma.
#    -x/--xalt not homozygous alternate
#    -n/--xref not homozygous reference   
#    -u/--unfl not filtered just output genotype data
#    -d/--dp    <optional>    The user may optionally specify a minimum depth of coverage for each individual
#    -g/--gq    <optional>    The user may optionally specify a minimum genotyping quality for each individual
#    -q/--mq    <optional>    The user may optionally specify a maximum for the ratio of MQ0/DP for each variant
#    -m/--maf    <optional>    The user may optionally specify a maximum alternate allele frequency for each variant - compared to both 1KG and ESP-GO
#
# It is not necessary to specify all of the samples in the vcf, if only a subset is supplied only that subset will be used for filtering and output to the tab delimited file; the output vcf will contain all samples. If no samples are specified then all variants will be returned
# If a variant has multiple alternate alleles then the script will, with some limitations, iterate across possible combinations of these for heterozygous and homozygous alternate. Limitations:
#            The homozygous alternate allele will be the same for all individuals (i.e. if two individual are 2/2 and 3/3, this will not be kept)
#            The heterozygous genotype will always contain the homozyogus alternate genotype allele (i.e. if the homozygous atlernate is 3/3, a variant with the heterozygote 1/2 will not be kept but one with 1/3 or 2/3 will)
# For example to run the autosomal dominant filter described above with the default quality filters:
#    ExmFilt.CustomFiltering.py -v MyData.vcf -o MyOutputFiles -e ProbandID,Parent1ID -r Parent2ID

import os
from optparse import OptionParser
parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="Input VCF file", metavar="VCFfile")
parser.add_option("-o", "--out", dest="OutputFileName",help="Base of names for output files", metavar="OutputFileName")

# Assign samples to heterozygous, homozygous or reference
parser.add_option("-a", "--alt", dest="homozygousgenotypes",help="Samples that should be homozygous for alternate allele, e.g. 1/1, 2/2, separated by commas", metavar="homozygousgenotypes")
parser.add_option("-e", "--het", dest="heterozygousgenotypes",help="Samples that should be heterozygous, e.g. 0/1, separated by commas", metavar="heterozygousgenotypes")
parser.add_option("-r", "--ref", dest="referencegenotypes",help="Samples that should be homozygous for reference allele, i.e. 0/0, separated by commas", metavar="referencegenotypes")
parser.add_option("-x", "--xalt", dest="notalternate", help="Samples that should not be homozygous for alternate allele, i.e. maybe 0/0 or 0/1, separated by commas", metavar="notalternate")
parser.add_option("-n", "--xref", dest="notreference",help="Samples that should not be homozygous for reference allele, i.e. maybe  0/1 or 1/1, separated by commas", metavar="notreference")
parser.add_option("-u", "--unfl", dest="notfiltered",help="Samples that should not be used for filtereing but output the genoytpe data", metavar="notfiltered")

# Additional Parameters to Adjust filtering thresholds
parser.set_defaults(MQ0threshold="0.05")
parser.set_defaults(DPthreshold="5")
parser.set_defaults(MAFthreshold="0.01")
parser.set_defaults(GQthreshold="40")


parser.add_option("-d", "--dp", dest="DPthreshold",help="minimum depth of coverage", metavar="DPthreshold")
parser.add_option("-m", "--mq", dest="MQ0threshold",help="maximum for MQ value, given as a decimal fraction (e.g. 0.05 = 5% of reads for a variant can have MQ0) - fraction of reads with quality 0", metavar="MQ0")
parser.add_option("-f", "--maf", dest="MAFthreshold",help="maximum MAF", metavar="MAFthreshold")
parser.add_option("-g", "--gq", dest="GQthreshold",help="minimum for GQ value", metavar="GQthreshold")

parser.add_option("-P", "--highpathogenicity", action='store_true', dest="highpatho", help="Filter using high pathogenicity predictions")
parser.add_option("-p", "--lowpathogenicity", action='store_true', dest="lowpatho", help="Filter using low pathogenicity predictions")
parser.add_option("-X", "--xlinked", action='store_true', dest="xlink", help="Only output X chromosome")


parser.add_option("-z", "--debug", action='store_true', dest="debug", help="Output a boolean decision file for debugging purposes")

(options, args) = parser.parse_args()


#Assign input and output files
VCF=open(options.VCFfile,'r')
BaseName=str(options.OutputFileName)
TabOutputFilename=BaseName+'.tsv'
VcfOutputFilename=BaseName+'.vcf'
LogOutputFilename=BaseName+'.log'
PassOutputFilename=BaseName+'.boolean.log'
Output=open(TabOutputFilename,'w')
Outvcf=open(VcfOutputFilename,'w')
Outlog=open(LogOutputFilename,'w')
HighPatho=options.highpatho
LowPatho=options.lowpatho
XLink=options.xlink
DeBug=options.debug
if DeBug:
    OutPass=open(PassOutputFilename,'w')

# Assign filter variables from optparse
MQ0filter=float(options.MQ0threshold)
MAFfilter=float(options.MAFthreshold)
DPfilter=int(options.DPthreshold)
GQfilter=int(options.GQthreshold)

# Define sample names using user-defined parameters
HomozygousSamples=options.homozygousgenotypes
HeterozygousSamples=options.heterozygousgenotypes
ReferenceSamples=options.referencegenotypes
NotAlternateSamples=options.notalternate
NotReferenceSamples=options.notreference
NotFilteredSamples=options.notfiltered
HomozygousSampleList=[]
HeterozygousSampleList=[]
ReferenceSampleList=[]
NotAlternateSampleList=[]
NotReferenceSampleList=[]
NotFilteredSampleList=[]


if HomozygousSamples is not None:
    HomozygousSampleList=HomozygousSamples.upper().split(',')
if HeterozygousSamples is not None:
    HeterozygousSampleList=HeterozygousSamples.upper().split(',')
if ReferenceSamples is not None:
    ReferenceSampleList=ReferenceSamples.upper().split(',')
if NotReferenceSamples is not None:
    NotReferenceSampleList=NotReferenceSamples.upper().split(',')
if NotAlternateSamples is not None:
    NotAlternateSampleList=NotAlternateSamples.upper().split(',')
if NotFilteredSamples is not None:
    NotFilteredSampleList=NotFilteredSamples.upper().split(',')

#start log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("De novo filtering log: "+TimeNow+"\n")
Outlog.write("VCF: "+str(options.VCFfile)+"\n")

Outlog.write("Genotype Filters: \n")
Outlog.write("Homozygous Samples: "+str(HomozygousSamples)+"\n")
Outlog.write("Heterozygous Samples: "+str(HeterozygousSamples)+"\n")
Outlog.write("Reference Samples: "+str(ReferenceSamples)+"\n")
Outlog.write("Not-Reference Samples: "+str(NotReferenceSamples)+"\n")
Outlog.write("Not-Alternate Samples: "+str(NotAlternateSamples)+"\n")
Outlog.write("Unfiltered Samples: "+str(NotFilteredSamples)+"\n")

Outlog.write("Variant Filters: \n")
Outlog.write("\t MQ0/DP maximum: "+str(MQ0filter)+"\n")
Outlog.write("\t 1000 genomes alternate allele frequency maximum: "+str(MAFfilter)+"\n")
Outlog.write("\t GO ESP alternate allele frequency maximum: "+str(MAFfilter)+"\n")

Outlog.write("Individual Sample Filters: \n")
Outlog.write("\t Minimum Depth of Coverage: "+str(DPfilter)+"\n")
Outlog.write("\t Minimum Genotyping Quality(GQ): "+str(GQfilter)+"\n")
if not LowPatho and not HighPatho:
    Outlog.write("\t Predicted Pathogenecity Filters: None \n")
if HighPatho:
    Outlog.write("\t Predicted Pathogenecity Filters: both SIFT=D AND PP2=D or P, or CADD >= 20 \n")
if LowPatho:
    Outlog.write("\t Predicted Pathogenecity Filters: SIFT=D or PP2=D or P or CADD >= 10 \n")

#Start output table
headerlist=['Chromosome','Position','ID','REF','ALT','Gene','VariantFunction','VariantClass','AAchange','AlleleFrequency.1KG','AlleleFrequency.ESP','SIFTprediction','PP2prediction','MTprediction','CADDscore']+HomozygousSampleList+HeterozygousSampleList+ReferenceSampleList+NotReferenceSampleList+NotAlternateSampleList+NotFilteredSampleList
Output.write("\t".join(headerlist)+"\n")

# Read VCF file
NameToColumn={}
ColumnToName={}
OrigCount=0
FiltCount=0

for line in VCF:
    # Map column name to number, and then find column numbers of each set of trios
    if '#' in line:
        Outvcf.write(line)
    if '#CHROM' in line:
        ColumnAndNumber=enumerate(line.strip().upper().split("\t"))
        for i in ColumnAndNumber:
            NameToColumn[i[1]]=i[0]
            ColumnToName[i[0]]=i[1]
        HomozygousColumnNumber=[ NameToColumn[i] for i in HomozygousSampleList ]
        HeterozygousColumnNumber=[ NameToColumn[i] for i in HeterozygousSampleList ]
        ReferenceColumnNumber=[ NameToColumn[i] for i in ReferenceSampleList ]
        NotReferenceColumnNumber=[ NameToColumn[i] for i in NotReferenceSampleList ]
        NotAlternateColumnNumber=[ NameToColumn[i] for i in NotAlternateSampleList ]
        NotFilteredColumnNumber=[ NameToColumn[i] for i in NotFilteredSampleList ]
    if '#' not in line:
        linelist=line.split("\t", 1)
        Chrom=linelist[0]
        PassX=True
        if Chrom!="X" and XLink:
            PassX=False
    if '#' not in line and PassX:
        # Variant must first pass 1KG and GO-ESP frequencies, MQ0 threshold, and be exonic
        OrigCount=OrigCount+1
        linelist=line.split("\t")
        SampleNum=len(linelist)-9
        
        REFlist=linelist[4]
        REFlist=REFlist.split(",")
        
        QUAL=linelist[5]
        INFOstring=linelist[7]
        INFOcolumnList=INFOstring.split(";")
        INFOdict={}
        for element in INFOcolumnList:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                INFOdict[FieldName]=FieldValue
        
        # Get variant data
        
        
        QDnumber=float(INFOdict.get('QD',0))
        DPnumber=float(INFOdict.get('DP',0))
        MQ0number=float(INFOdict.get('MQ0',0))
        GeneName=INFOdict.get('GeneName','NA')
        VariantFunction=INFOdict.get('VarFunc','none')
        
        KGscoreList=str(INFOdict.get('1KGfreq',0))
        KGscoreList=KGscoreList.split(',')
        ESPscoreList=str(INFOdict.get('ESPfreq',0))
        ESPscoreList=ESPscoreList.split(',')
        
        VariantClassList=INFOdict.get('VarClass','none')
        VariantClassList=VariantClassList.split(',')
        AAchangeList=INFOdict.get('AAChange','NA')
        AAchangeList=AAchangeList.split(',')
        SIFTscoreList=str(INFOdict.get('SIFTscr',0))
        SIFTscoreList=SIFTscoreList.split(',')
        SIFTpredictionList=INFOdict.get('SIFTprd','NA')
        SIFTpredictionList=SIFTpredictionList.split(',')
        PP2scoreList=str(INFOdict.get('PP2.hvar.scr',0))
        PP2scoreList=PP2scoreList.split(',')
        PP2predictionList=INFOdict.get('PP2.hvar.prd','NA')
        PP2predictionList=PP2predictionList.split(',')
        MTscoreList=str(INFOdict.get('MutTscr',0))
        MTscoreList=MTscoreList.split(',')
        MTpredictionList=INFOdict.get('MutTprd','NA')
        MTpredictionList=MTpredictionList.split(',')
        GERPscoreList=str(INFOdict.get('GERP',0))
        GERPscoreList=GERPscoreList.split(',')
        PHYLOPscoreList=str(INFOdict.get('PhyloP',0))
        PHYLOPscoreList=PHYLOPscoreList.split(',')
        CADDscoreList=str(INFOdict.get('CADDphred',0))
        CADDscoreList=CADDscoreList.split(',')
        
        # Check if MQ0 passes threshold
        PassMQ=False
        MQ0Fraction=MQ0number/DPnumber
        if MQ0Fraction <= MQ0filter:
            PassMQ=True
         
        # Get number of alternate alleles
        AltAlls=linelist[4]
        AltAlls=AltAlls.split(",")
        AltNum=len(AltAlls)
        
        # Get Sample Qualities
        HomozygousQualityString=[ linelist[i].strip() for i in HomozygousColumnNumber ]
        HeterozygousQualityString=[ linelist[i].strip() for i in HeterozygousColumnNumber ]
        ReferenceQualityString=[ linelist[i].strip() for i in ReferenceColumnNumber ]
        NotReferenceQualityString=[ linelist[i].strip() for i in NotReferenceColumnNumber ]
        NotAlternateQualityString=[ linelist[i].strip() for i in NotAlternateColumnNumber ]
        NotFilteredQualityString=[ linelist[i].strip() for i in NotFilteredColumnNumber ]
        
        # Split Individual Quality strings for each sample
        HomozygousQualityList=[ i.split(':') for i in HomozygousQualityString ]
        HeterozygousQualityList=[ i.split(':') for i in HeterozygousQualityString ]
        ReferenceQualityList=[ i.split(':') for i in ReferenceQualityString ]
        NotReferenceQualityList=[ i.split(':') for i in NotReferenceQualityString ]
        NotAlternateQualityList=[ i.split(':') for i in NotAlternateQualityString ]
        NotFilteredQualityList=[ i.split(':') for i in NotFilteredQualityString ]
        
        # Define GT
        HomozygousGT=[ HomozygousQualityList[i][0] for i in range(0,len(HomozygousQualityString)) ]
        HeterozygousGT=[ HeterozygousQualityList[i][0] for i in range(0,len(HeterozygousQualityString)) ]
        ReferenceGT=[ ReferenceQualityList[i][0] for i in range(0,len(ReferenceQualityString)) ]
        NotReferenceGT=[ NotReferenceQualityList[i][0] for i in range(0,len(NotReferenceQualityString)) ]
        NotAlternateGT=[ NotAlternateQualityList[i][0] for i in range(0,len(NotAlternateQualityString)) ]
        NotFilteredGT=[ NotFilteredQualityList[i][0] for i in range(0,len(NotFilteredQualityString)) ]
        
        # First check if all genotypes are present
        
        
        PassFunction=False
        if VariantFunction=='exonic' or VariantFunction=='splicing' or VariantFunction=="exonic,splicing" or VariantFunction=='none':
            PassFunction=True
        
        # Check if Variant class passes
        if './.' not in HomozygousGT and './.' not in HeterozygousGT and './.' not in ReferenceGT and './.' not in NotReferenceGT and './.' not in NotAlternateGT and PassFunction:
            
            # Define DP
            HomozygousDP=[ str(HomozygousQualityList[i][2]) for i in range(0,len(HomozygousQualityString)) ]
            for i in range(0,len(HomozygousDP)):
                if HomozygousDP[i]==".":
                    HomozygousDP[i]="0"
                HomozygousDP[i]=float(HomozygousDP[i])
            HeterozygousDP=[ str(HeterozygousQualityList[i][2]) for i in range(0,len(HeterozygousQualityString)) ]
            for i in range(0,len(HeterozygousDP)):
                if HeterozygousDP[i]==".":
                    HeterozygousDP[i]="0"
                HeterozygousDP[i]=float(HeterozygousDP[i])
            ReferenceDP=[ str(ReferenceQualityList[i][2]) for i in range(0,len(ReferenceQualityString)) ]
            for i in range(0,len(ReferenceDP)):
                if ReferenceDP[i]==".":
                    ReferenceDP[i]="0"
                ReferenceDP[i]=float(ReferenceDP[i])
            NotReferenceDP=[ float(NotReferenceQualityList[i][2]) for i in range(0,len(NotReferenceQualityString)) ]
            NotAlternateDP=[ float(NotAlternateQualityList[i][2]) for i in range(0,len(NotAlternateQualityString)) ]
            # Define GQ
            HomozygousGQ=[ float(HomozygousQualityList[i][3]) for i in range(0,len(HomozygousQualityString)) ]
            HeterozygousGQ=[ float(HeterozygousQualityList[i][3]) for i in range(0,len(HeterozygousQualityString)) ]
            ReferenceGQ=[ float(ReferenceQualityList[i][3]) for i in range(0,len(ReferenceQualityString)) ]
            NotReferenceGQ=[ float(NotReferenceQualityList[i][3]) for i in range(0,len(NotReferenceQualityString)) ]
            NotAlternateGQ=[ float(NotAlternateQualityList[i][3]) for i in range(0,len(NotAlternateQualityString)) ]
            
            # Check to see if depth passes
            PassDP=False
            if all(i >= DPfilter for i in HomozygousDP) and all(i >= DPfilter for i in HeterozygousDP) and all(i >= DPfilter for i in ReferenceDP):
                PassDP=True
            # Check to see if genotype quality passes
            PassGQ=False
            if all(i >= GQfilter for i in HomozygousGQ) and all(i >= GQfilter for i in HeterozygousGQ) and all(i >= GQfilter for i in ReferenceGQ):
                PassGQ=True
            # Check to see if genotypes pass, iterate across multiple alt alleles if necessary
            AltRng=range(0, AltNum)
            for altnum in AltRng:
                #check variant class
                PassClass=False
                cltnum=min(len(VariantClassList)-1, altnum)
                VariantClass=str(VariantClassList[cltnum])
                if VariantClass != 'synonymousSNV':
                    PassClass=True
                REF=str(REFlist[altnum])
                
                # Check if KG passes threshold
                PassKG=False
                cltnum=min(len(KGscoreList)-1, altnum)
                KGscore=str(KGscoreList[cltnum])
                if KGscore == ".":
                    KGscore="0"
                KGscore=float(KGscore)
                if KGscore <= MAFfilter:
                    PassKG=True
                
                # Check if ESP passes threshold
                PassESP=False
                cltnum=min(len(ESPscoreList)-1, altnum)
                ESPscore=str(ESPscoreList[cltnum])
                if ESPscore == ".":
                    ESPscore="0"
                ESPscore=float(ESPscore)
                if ESPscore <= MAFfilter:
                    PassESP=True
                
                
                #check GT
                PassGT=False
                AltAll=str(altnum+1)
                if all(i.count(AltAll)==2 for i in HomozygousGT) and all(i.count(AltAll)==1 for i in HeterozygousGT) and all(i.count(AltAll)==0 for i in ReferenceGT) and all(i.count(AltAll)>0 for i in NotReferenceGT) and all(i.count(AltAll)<2 for i in NotAlternateGT):
                    PassGT=True
                
                PassPatho=False
                cltnum=min(len(SIFTpredictionList)-1, altnum)
                SIFTprediction=SIFTpredictionList[cltnum]
                cltnum=min(len(PP2predictionList)-1, altnum)
                PP2prediction=PP2predictionList[cltnum]
                cltnum=min(len(CADDscoreList)-1, altnum)
                CADDscore=str(CADDscoreList[cltnum])
                if CADDscore == ".":
                    CADDscore="0"
                CADDscore=float(CADDscore)
                if HighPatho and (SIFTprediction=="D" and PP2prediction!="B"):
                    PassPatho=True
                if HighPatho and CADDscore>=20:
                    PassPatho=True
                if LowPatho and (SIFTprediction=="D" or PP2prediction=="D" or PP2prediction=="P" or CADDscore>=15):
                    PassPatho=True
                if not LowPatho and not HighPatho:
                    PassPatho=True
                
                
                cltnum=min(len(ESPscoreList)-1, altnum)
                ESPscore=str(ESPscoreList[cltnum])
                cltnum=min(len(KGscoreList)-1, altnum)
                KGscore=str(KGscoreList[cltnum])
                cltnum=min(len(CADDscoreList)-1, altnum)
                CADDscore=str(CADDscoreList[cltnum])
                cltnum=min(len(AAchangeList)-1, altnum)
                AAchange=AAchangeList[cltnum]
                cltnum=min(len(MTpredictionList)-1, altnum)
                MTprediction=MTpredictionList[cltnum]
                cltnum=min(len(PHYLOPscoreList)-1, altnum)
                PHYLOPscore=PHYLOPscoreList[cltnum]
                #output
                if DeBug:
                    OutPass.write("\t"+linelist[0]+" "+linelist[1]+" "+linelist[3]+" "+linelist[4]+" "+str(altnum)+" "+str(PassKG)+" "+str(PassESP)+" "+str(PassMQ)+" "+str(PassFunction)+" "+str(PassClass)+" "+str(PassGT)+" "+str(PassDP)+" "+str(PassGQ)+" "+str(PassPatho)+"\n")
                # If all pass then output line
                if PassKG and PassESP and PassMQ and PassGT and PassDP and PassGQ and PassPatho:
                    OutputList=linelist[0:4]+[REF,GeneName,VariantFunction,VariantClass,AAchange,KGscore,ESPscore,SIFTprediction,PP2prediction,MTprediction,CADDscore]+HomozygousGT+HeterozygousGT+ReferenceGT+NotReferenceGT+NotAlternateGT+NotFilteredGT
                    OutputList= [ str(i) for i in OutputList ]
                    OutputString="\t".join(OutputList)
                    Output.write(OutputString+"\n")
                    Outvcf.write(line)
                    FiltCount=FiltCount+1
Outlog.write("Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("Number of variants selected: "+str(FiltCount)+"\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Filtering Finished: "+TimeNow+"\n")
print "Filtering finished "+TimeNow