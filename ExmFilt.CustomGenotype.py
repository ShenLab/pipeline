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
#    -m/--mq    <optional>    The user may optionally specify a maximum for the ratio of MQ0/DP for each variant
#    -f/--maf    <optional>    The user may optionally specify a maximum alternate allele frequency for each variant - compared to both 1KG and ESP-GO
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
parser.set_defaults(DPthreshold="3")
parser.set_defaults(MAFthreshold="0.01")
parser.set_defaults(GQthreshold="30")
parser.set_defaults(CHTthreshold="0.2")
parser.set_defaults(HetAllCountthreshold="3")
parser.set_defaults(HomAllFracthreshold="0.8")

## b        k l    p q  s t   w  y z 

parser.add_option("-d", "--dp", dest="DPthreshold",help="minimum depth of coverage", metavar="DPthreshold")
parser.add_option("-m", "--mq", dest="MQ0threshold",help="maximum for MQ value, given as a decimal fraction (e.g. 0.05 = 5% of reads for a variant can have MQ0) - fraction of reads with quality 0", metavar="MQ0")
parser.add_option("-f", "--maf", dest="MAFthreshold",help="maximum MAF", metavar="MAFthreshold")
parser.add_option("-g", "--gq", dest="GQthreshold",help="minimum for GQ value", metavar="GQthreshold")
parser.add_option("-c", "--cht", dest="CHTthreshold",help="maximum frequency of allele in cohort of samples in the vcf", metavar="CHTthreshold")
parser.add_option("-i", "--aac", dest="HetAllCountthreshold",help="minimum alternate allele count in heterozygous calls", metavar="HetAllCountthreshold")
parser.add_option("-j", "--haf", dest="HomAllFracthreshold",help="Fraction of allele counts for homozygous allele", metavar="HomAllFracthreshold")

parser.add_option("-P", "--nopathogenicity", action='store_true', dest="nopatho", help="Do not filter using pathogenicity predictions")

parser.add_option("-X", "--xlinked", action='store_true', dest="xlink", help="Only output X chromosome")
parser.add_option("-D", "--denovo", action='store_true', dest="denovo", help="Use higher stringency filters for de novo")
parser.add_option("-Q", "--noqual", action='store_true', dest="noqual", help="Do NOT filter by vcf FILTER field")


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
NoPatho=options.nopatho
XLink=options.xlink
DeNovo=options.denovo
DeBug=options.debug
NoQualityFilter=options.noqual
if DeBug:
    OutPass=open(PassOutputFilename,'w')

# Assign filter variables from optparse
HetAllCountFilter=int(options.HetAllCountthreshold)
HomAllFrac=float(options.HomAllFracthreshold)
MQ0filter=float(options.MQ0threshold)
MAFfilter=float(options.MAFthreshold)
DPfilter=int(options.DPthreshold)
GQfilter=int(options.GQthreshold)
CHTfilter=float(options.CHTthreshold)
if DeNovo:
    HetAllCountFilter=6
    HomAllFrac=0.98
    MAFfilter=0.001
    CHTfilter=0.01
    DPfilter=5

# Define sample names using user-defined parameters
AlternateSamples=options.homozygousgenotypes
HeterozygousSamples=options.heterozygousgenotypes
ReferenceSamples=options.referencegenotypes
NotAlternateSamples=options.notalternate
NotReferenceSamples=options.notreference
NotFilteredSamples=options.notfiltered
AlternateSampleList=[]
HeterozygousSampleList=[]
ReferenceSampleList=[]
NotAlternateSampleList=[]
NotReferenceSampleList=[]
NotFilteredSampleList=[]


if AlternateSamples is not None:
    AlternateSampleList=AlternateSamples.upper().split(',')
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
Outlog.write("Filtering log: "+TimeNow+"\n")
Outlog.write("VCF: "+str(options.VCFfile)+"\n")
Outlog.write("OutputName: "+str(options.OutputFileName)+"\n")
Outlog.write("\n")
Outlog.write("Genotype Filters: \n")
Outlog.write("\t Homozygous Alternate Samples: "+str(AlternateSamples)+"\n")
Outlog.write("\t Heterozygous Samples: "+str(HeterozygousSamples)+"\n")
Outlog.write("\t Homozygous Reference Samples: "+str(ReferenceSamples)+"\n")
Outlog.write("\t Not-Reference Samples: "+str(NotReferenceSamples)+"\n")
Outlog.write("\t Not-Alternate Samples: "+str(NotAlternateSamples)+"\n")
Outlog.write("\t Unfiltered Samples: "+str(NotFilteredSamples)+"\n")
Outlog.write("\n")
Outlog.write("Variant Filters: \n")
Outlog.write("\t MQ0/DP maximum: "+str(MQ0filter)+"\n")
Outlog.write("\t 1000 genomes alternate allele frequency maximum: "+str(MAFfilter)+"\n")
Outlog.write("\t GO ESP alternate allele frequency maximum: "+str(MAFfilter)+"\n")
Outlog.write("\t Within VCF allele frequency maximum: "+str(CHTfilter)+"\n")
Outlog.write("\n")
Outlog.write("Individual Sample Filters: \n")
Outlog.write("\t Minimum Depth of Coverage: "+str(DPfilter)+"\n")
Outlog.write("\t Minimum Genotyping Quality(GQ): "+str(GQfilter)+"\n")
Outlog.write("\t Minimum Heterozygous allele count: "+str(HetAllCountFilter)+"\n")
Outlog.write("\t Minimum Homozygous allele fraction: "+str(HomAllFrac)+"\n")
if NoPatho:
    Outlog.write("\t No Predicted Pathogenecity Filters")
if not NoPatho:
    Outlog.write("\t Predicted Pathogenecity Filters: at least SIFT=D or PP2=D/P or CADD >= 15; all InDels and Nonsense \n")
Outlog.write("\n")

#Start output table
headerlist=['Chromosome','Position','ID','REF','ALT','Gene','VariantFunction','VariantClass','AAchange','AlleleFrequency.1KG','AlleleFrequency.ESP','SIFTprediction','PP2prediction','MTprediction','GERP++','CADDscore','PredictionSummary']+AlternateSampleList+HeterozygousSampleList+ReferenceSampleList+NotReferenceSampleList+NotAlternateSampleList+NotFilteredSampleList+['AlternateAlleles', 'FILTER', 'INFO']+AlternateSampleList+HeterozygousSampleList+ReferenceSampleList+NotReferenceSampleList+NotAlternateSampleList+NotFilteredSampleList
Output.write("\t".join(headerlist)+"\n")

# Read VCF file
NameToColumn={}
ColumnToName={}
OrigCount=0
FiltCount=[0,0]
MissenseCount=[0,0]
NonFrameshiftCount=[0,0]
FrameshiftCount=[0,0]
NonsenseCount=[0,0]
SplicingCount=[0,0]
UnknownCount=[0,0]

MissenseClass=['nonsynonymousSNV','unknown']
NonFrameshiftClass=['nonframeshiftdeletion','nonframeshiftinsertion','nonframeshiftsubstitution']
FrameshiftClass=['frameshiftdeletion','frameshiftinsertion','frameshiftsubstitution']
InDelClass=NonFrameshiftClass+FrameshiftClass
NonsenseClass=['stopgain','stoploss']
SplicingClass=['splicing','exonic,splicing']

#BadSnpFilters=['QD_Bad_SNP','FS_Bad_SNP','FS_Mid_SNP;QD_Mid_SNP','SnpCluster','StandardFilters']
#BadInDFilters=['FSBias_Indel','LowQD_Indel','RPBias_Indel','SnpCluster','StandardFilters']
BadSnpFilters=['QD_Bad_SNP','FS_Bad_SNP','FS_Mid_SNP;QD_Mid_SNP','StandardFilters']
BadInDFilters=['FSBias_Indel','LowQD_Indel','RPBias_Indel','StandardFilters']



for line in VCF:
    # Map column name to number, and then find column numbers of each set of trios
    if '#' in line:
        Outvcf.write(line)
    if '#CHROM' in line:
        ColumnAndNumber=enumerate(line.strip().upper().split("\t"))
        for i in ColumnAndNumber:
            NameToColumn[i[1]]=i[0]
            ColumnToName[i[0]]=i[1]
        AlternateColumnNumber=[ NameToColumn[i] for i in AlternateSampleList ]
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
        
        REFstring=linelist[4]
        REFlist=REFstring.split(",")
        
        VariantFilter=linelist[6]
        VariantFilterList=VariantFilter.split(";")
        INFOstring=linelist[7]
        INFOcolumnList=INFOstring.split(";")
        INFOdict={}
        for element in INFOcolumnList:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                INFOdict[FieldName]=FieldValue
        
        # Get variant data
        
        DPnumber=float(INFOdict.get('DP','.'))
        MQ0number=float(INFOdict.get('MQ0','.'))
        GeneName=INFOdict.get('GeneName','.')
        VariantFunction=INFOdict.get('VarFunc','none')
        
        KGFreqList=str(INFOdict.get('1KGfreq','.'))
        KGFreqList=KGFreqList.split(',')
        ESPFreqList=str(INFOdict.get('ESPfreq','.'))
        ESPFreqList=ESPFreqList.split(',')
        VCFFreqList=str(INFOdict.get('AF',0))
        VCFFreqList=VCFFreqList.split(",")
        
        VariantClassList=INFOdict.get('VarClass','none')
        VariantClassList=VariantClassList.split(',')
        AAchangeList=INFOdict.get('AAChange','.')
        AAchangeList=AAchangeList.split(',')
        SIFTscoreList=str(INFOdict.get('SIFTscr','.'))
        SIFTscoreList=SIFTscoreList.split(',')
        SIFTpredictionList=INFOdict.get('SIFTprd','.')
        SIFTpredictionList=SIFTpredictionList.split(',')
        PP2scoreList=str(INFOdict.get('PP2.hvar.scr','.'))
        PP2scoreList=PP2scoreList.split(',')
        PP2predictionList=INFOdict.get('PP2.hvar.prd','.')
        PP2predictionList=PP2predictionList.split(',')
        MTscoreList=str(INFOdict.get('MutTscr','.'))
        MTscoreList=MTscoreList.split(',')
        MTpredictionList=INFOdict.get('MutTprd','.')
        MTpredictionList=MTpredictionList.split(',')
        GERPscoreList=str(INFOdict.get('GERP','.'))
        GERPscoreList=GERPscoreList.split(',')
        CADDscoreList=str(INFOdict.get('CADDphred','.'))
        CADDscoreList=CADDscoreList.split(',')
        
        #Variant level checks
        
        # Get number of alternate alleles
        AltAlls=linelist[4]
        AltAlls=AltAlls.split(",")
        AltNum=len(AltAlls)
        
        # Get Sample Qualities
        AlternateQualityString=[ linelist[i].strip() for i in AlternateColumnNumber ]
        HeterozygousQualityString=[ linelist[i].strip() for i in HeterozygousColumnNumber ]
        ReferenceQualityString=[ linelist[i].strip() for i in ReferenceColumnNumber ]
        NotReferenceQualityString=[ linelist[i].strip() for i in NotReferenceColumnNumber ]
        NotAlternateQualityString=[ linelist[i].strip() for i in NotAlternateColumnNumber ]
        NotFilteredQualityString=[ linelist[i].strip() for i in NotFilteredColumnNumber ]
        
        # Split Individual Quality strings for each sample
        AlternateQualityList=[ i.split(':') for i in AlternateQualityString ]
        HeterozygousQualityList=[ i.split(':') for i in HeterozygousQualityString ]
        ReferenceQualityList=[ i.split(':') for i in ReferenceQualityString ]
        NotReferenceQualityList=[ i.split(':') for i in NotReferenceQualityString ]
        NotAlternateQualityList=[ i.split(':') for i in NotAlternateQualityString ]
        NotFilteredQualityList=[ i.split(':') for i in NotFilteredQualityString ]
        
        # Define GT
        AlternateGT=[ AlternateQualityList[i][0] for i in range(0,len(AlternateQualityString)) ]
        HeterozygousGT=[ HeterozygousQualityList[i][0] for i in range(0,len(HeterozygousQualityString)) ]
        ReferenceGT=[ ReferenceQualityList[i][0] for i in range(0,len(ReferenceQualityString)) ]
        NotReferenceGT=[ NotReferenceQualityList[i][0] for i in range(0,len(NotReferenceQualityString)) ]
        NotAlternateGT=[ NotAlternateQualityList[i][0] for i in range(0,len(NotAlternateQualityString)) ]
        NotFilteredGT=[ NotFilteredQualityList[i][0] for i in range(0,len(NotFilteredQualityString)) ]
        
        # Check if MQ0 passes threshold
        PassMQ=False
        MQ0Fraction=MQ0number/DPnumber
        if MQ0Fraction <= MQ0filter:
            PassMQ=True
         
        # Check if all genotypes are present
        PassGeno=False
        if './.' not in AlternateGT and './.' not in HeterozygousGT and './.' not in ReferenceGT and './.' not in NotReferenceGT and './.' not in NotAlternateGT:
            PassGeno=True
        
        # Check if Variant class passes
        PassFunction=False
        if VariantFunction=='exonic' or VariantFunction in SplicingClass or VariantFunction=='none':
            PassFunction=True
        
        if PassGeno and PassFunction and PassMQ :
            
            # Define DP
            AlternateDP=[ str(AlternateQualityList[i][2]) for i in range(0,len(AlternateQualityString)) ]
            for i in range(0,len(AlternateDP)):
                if AlternateDP[i]==".":
                    AlternateDP[i]="0"
                AlternateDP[i]=float(AlternateDP[i])
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
            for i in range(0,len(NotReferenceDP)):
                if NotReferenceDP[i]==".":
                    NotReferenceDP[i]="0"
                NotReferenceDP[i]=float(NotReferenceDP[i])
            NotAlternateDP=[ float(NotAlternateQualityList[i][2]) for i in range(0,len(NotAlternateQualityString)) ]
            for i in range(0,len(NotAlternateDP)):
                if NotAlternateDP[i]==".":
                    NotAlternateDP[i]="0"
                NotAlternateDP[i]=float(NotAlternateDP[i])
            
            # Define GQ
            AlternateGQ=[ float(AlternateQualityList[i][3]) for i in range(0,len(AlternateQualityString)) ]
            HeterozygousGQ=[ float(HeterozygousQualityList[i][3]) for i in range(0,len(HeterozygousQualityString)) ]
            ReferenceGQ=[ float(ReferenceQualityList[i][3]) for i in range(0,len(ReferenceQualityString)) ]
            NotReferenceGQ=[ float(NotReferenceQualityList[i][3]) for i in range(0,len(NotReferenceQualityString)) ]
            NotAlternateGQ=[ float(NotAlternateQualityList[i][3]) for i in range(0,len(NotAlternateQualityString)) ]
            
            # Check to see if depth passes
            PassDP=False
            if all(i >= DPfilter for i in AlternateDP) and all(i >= DPfilter for i in HeterozygousDP) and all(i >= DPfilter for i in ReferenceDP):
                PassDP=True
            # Check to see if genotype quality passes
            PassGQ=False
            if all(i >= GQfilter for i in AlternateGQ) and all(i >= GQfilter for i in HeterozygousGQ) and all(i >= GQfilter for i in ReferenceGQ):
                PassGQ=True
            # Check to see if genotypes pass, iterate across multiple alt alleles if necessary
            AltRng=range(0, AltNum)
            for altnum in AltRng:
                
                
                # Define Allele Count
                
                AlternateAC=[ AlternateQualityList[i][1] for i in range(0,len(AlternateQualityString)) ]
                AlternateAC=[ i.split(',') for i in AlternateAC ]
                AlternateTotal=[ sum(int(j) for j in i)  for i in AlternateAC ]
                AlternateAC=[ int(i[altnum+1]) for i in AlternateAC ]
                AlternateAAF=[ float(i)/float(max(j,1)) for i,j in zip(AlternateAC,AlternateTotal) ] #if Total is 0 then python throws and error, so set to 0, this will still leave the AAF as 0 as if the total is 0 then the AC must also be 0
                
                HeterozygousAC=[ HeterozygousQualityList[i][1] for i in range(0,len(HeterozygousQualityString)) ]
                HeterozygousAC=[ i.split(',') for i in HeterozygousAC ]
                HeterozygousAC=[ int(i[altnum+1]) for i in HeterozygousAC ]
                
                ReferenceAC=[ ReferenceQualityList[i][1] for i in range(0,len(ReferenceQualityString)) ]
                ReferenceAC=[ i.split(',') for i in ReferenceAC ]
                ReferenceTotal=[ sum(int(j) for j in i)  for i in ReferenceAC ]
                ReferenceAC=[ int(i[0]) for i in ReferenceAC ]
                ReferenceAAF=[ float(i)/float(max(j,1)) for i,j in zip(ReferenceAC,ReferenceTotal) ]
                
                NotReferenceAC=[ NotReferenceQualityList[i][1] for i in range(0,len(NotReferenceQualityString)) ]
                NotReferenceAC=[ i.split(',') for i in NotReferenceAC ]
                NotReferenceTotal=[ sum(int(j) for j in i)  for i in NotReferenceAC ]
                NotReferenceAC=[ int(i[altnum+1]) for i in NotReferenceAC ]
                NotReferenceAAF=[ float(i)/float(max(j,1)) for i,j in zip(NotReferenceAC,NotReferenceTotal) ]
                
                NotAlternateAC=[ NotAlternateQualityList[i][1] for i in range(0,len(NotAlternateQualityString)) ]
                NotAlternateAC=[ i.split(',') for i in NotAlternateAC ]
                NotAlternateTotal=[ sum(int(j) for j in i)  for i in NotAlternateAC ]
                NotAlternateAAC=[ int(i[altnum+1]) for i in NotAlternateAC ]
                NotAlternateRAC=[ int(i[0]) for i in NotAlternateAC ]
                NotAlternateRAF=[ float(i)/float(max(j,1)) for i,j in zip(NotAlternateRAC,NotAlternateTotal) ]
                
                #check variant class
                PassClass=False
                cltnum=min(len(VariantClassList)-1, altnum)
                VariantClass=str(VariantClassList[cltnum])
                if VariantClass != 'synonymousSNV':
                    PassClass=True
                
                # Check if KG passes threshold
                PassKG=False
                cltnum=min(len(KGFreqList)-1, altnum)
                KGFreq=str(KGFreqList[cltnum])
                if KGFreq == ".":
                    KGFreqtest=float(0)
                else:
                    KGFreqtest=float(KGFreq)
                if KGFreqtest <= MAFfilter:
                    PassKG=True
                
                # Check if ESP passes threshold
                PassESP=False
                cltnum=min(len(ESPFreqList)-1, altnum)
                ESPFreq=str(ESPFreqList[cltnum])
                if ESPFreq == ".":
                    ESPFreqtest=float(0)
                else:
                    ESPFreqtest=float(ESPFreq)
                if ESPFreqtest <= MAFfilter:
                    PassESP=True
                
                # Check if VCF passes threshold
                PassVCF=False
                cltnum=min(len(VCFFreqList)-1, altnum)
                VCFFreqtest=float(VCFFreqList[cltnum])
                if VCFFreqtest <= CHTfilter:
                    PassVCF=True
                
                #check GT
                PassGT=False
                AltAll=str(altnum+1)
                if all(i.count(AltAll)==2 for i in AlternateGT) and all(i.count(AltAll)==1 for i in HeterozygousGT) and all(i.count(AltAll)==0 for i in ReferenceGT) and all(i.count(AltAll)>0 for i in NotReferenceGT) and all(i.count(AltAll)<2 for i in NotAlternateGT):
                    PassGT=True
                
                #Check alternate allele counts and fractions
                passBasicACF=False
                if all( i >= HetAllCountFilter for i in HeterozygousAC) and all( i >= HomAllFrac for i in AlternateAAF) and all( i >= HomAllFrac for i in ReferenceAAF):
                    passBasicACF=True
                
                passNotRefAFC=False
                for NotRefGT in NotReferenceGT:
                    if NotRefGT.count(AltAll) == 1 and NotReferenceAC >= HetAllCountFilter:
                        passNotRefAFC=True
                    if NotRefGT.count(AltAll) == 2 and NotReferenceAAF >= AlternateAAF:
                        passNotRefAFC=True
                if not NotReferenceGT:
                    passNotRefAFC=True
                
                passNotAltAFC=False
                for NotAltGT in NotAlternateGT:
                    if NotAltGT.count(AltAll) == 1 and NotAlternateAAC >= HetAllCountFilter:
                        passNotAltAFC=True
                    if NotAltGT.count(AltAll) == 0 and NotAlternateRAF >= AlternateAAF:
                        passNotAltAFC=True
                if not NotAlternateGT:
                    passNotAltAFC=True
                
                PassAFC=False
                if passBasicACF and passNotRefAFC and passNotAltAFC:
                    PassAFC=True
                
                PassPatho=False
                cltnum=min(len(SIFTpredictionList)-1, altnum)
                SIFTprediction=SIFTpredictionList[cltnum]
                cltnum=min(len(PP2predictionList)-1, altnum)
                PP2prediction=PP2predictionList[cltnum]
                cltnum=min(len(CADDscoreList)-1, altnum)
                CADDscore=str(CADDscoreList[cltnum])
                if CADDscore == ".":
                    CADDscoretest=float(0)
                else:
                    CADDscoretest=float(CADDscore)
                
                PathoLevel="Low"
                if VariantClass in MissenseClass and NoPatho:
                    PassPatho=True
                if VariantClass in MissenseClass and (SIFTprediction=="D" or PP2prediction=="D" or PP2prediction=="P" or CADDscoretest>=15):
                    PathoLevel="Med"
                    PassPatho=True
                if VariantClass in MissenseClass and (SIFTprediction=="D" and PP2prediction=="D"):
                    PathoLevel="High"
                    PassPatho=True
                if VariantClass in MissenseClass and (CADDscoretest>=25):
                    PathoLevel="High"
                    PassPatho=True
                if VariantFunction in SplicingClass:
                    PathoLevel="High"
                    PassPatho=True
                if VariantClass in InDelClass or VariantClass in NonsenseClass:
                    PathoLevel="High"
                    PassPatho=True
                
                
                PassFilter=False
                if VariantClass in InDelClass and all( str(i) not in BadInDFilters for i in VariantFilterList):
                    PassFilter=True
                if VariantClass not in InDelClass and all( str(i) not in BadSnpFilters for i in VariantFilterList):
                    PassFilter=True
                if NoQualityFilter:
                    PassFilter=True

                
                cltnum=min(len(AAchangeList)-1, altnum)
                AAchange=AAchangeList[cltnum]
                cltnum=min(len(MTpredictionList)-1, altnum)
                MTprediction=MTpredictionList[cltnum]
                #output
                if DeBug:
                    OutPass.write("\t"+linelist[0]+" "+linelist[1]+" "+linelist[3]+" "+linelist[4]+" "+str(altnum)+" "+str(PassKG)+" "+str(PassESP)+" "+str(PassVCF)+" "+str(PassMQ)+" "+str(PassFunction)+" "+str(PassClass)+" "+str(PassGT)+" "+str(PassDP)+" "+str(PassGQ)+" "+str(PassPatho)+" "+str(PassFilter)+" "+str(PassAFC)+"\n")
                # If all pass then output line
                if PassKG and PassESP and PassVCF and PassGT and PassDP and PassGQ and PassPatho and PassFilter and PassAFC:
                    REF=str(REFlist[altnum])
                    cltnum=min(len(GERPscoreList)-1, altnum)
                    GERPscore=str(GERPscoreList[cltnum])
                    OutputList=linelist[0:4]+[REF,GeneName,VariantFunction,VariantClass,AAchange,KGFreq,ESPFreq,SIFTprediction,PP2prediction,MTprediction,GERPscore,CADDscore,PathoLevel]+AlternateGT+HeterozygousGT+ReferenceGT+NotReferenceGT+NotAlternateGT+NotFilteredGT+[REFstring,VariantFilter,INFOstring]+AlternateQualityString+HeterozygousQualityString+ReferenceQualityString+NotReferenceQualityString+NotAlternateQualityString+NotFilteredQualityString
                    OutputList= [ str(i) for i in OutputList ]
                    OutputString="\t".join(OutputList)
                    Output.write(OutputString+"\n")
                    Outvcf.write(line)
                    cntnum=0
                    if PathoLevel=="High":
                        cntnum=1
                    FiltCount[cntnum] = FiltCount[cntnum] + 1
                    if VariantClass in MissenseClass:
                        MissenseCount[cntnum] = MissenseCount[cntnum] + 1
                    if VariantClass in NonFrameshiftClass:
                        NonFrameshiftCount[cntnum] = NonFrameshiftCount[cntnum] + 1
                    if VariantClass in FrameshiftClass:
                        FrameshiftCount[cntnum] = FrameshiftCount[cntnum] + 1
                    if VariantClass in NonsenseClass:
                        NonsenseCount[cntnum] = NonsenseCount[cntnum] + 1
                    if VariantFunction in SplicingClass:
                        SplicingCount[cntnum] = SplicingCount[cntnum] + 1
                    if VariantClass=="Unknown":
                        UnknownCount[cntnum] = UnknownCount[cntnum] + 1
Outlog.write("---------------------------------------------------------------------\n")
Outlog.write("Results: \n")
Outlog.write("\t Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("\t Number of variants selected (Med/High/Total): "+str(FiltCount[0])+"/"+str(FiltCount[1])+"/"+str(sum(FiltCount))+"\n")
Outlog.write("\t Number of Missense variants selected (Med/High/Total): "+str(MissenseCount[0])+"/"+str(MissenseCount[1])+"/"+str(sum(MissenseCount))+"\n")
Outlog.write("\t Number of Nonsense variants selected (Med/High/Total): "+str(NonsenseCount[0])+"/"+str(NonsenseCount[1])+"/"+str(sum(NonsenseCount))+"\n")
Outlog.write("\t Number of NonFrameshift InDels selected (Med/High/Total): "+str(NonFrameshiftCount[0])+"/"+str(NonFrameshiftCount[1])+"/"+str(sum(NonFrameshiftCount))+"\n")
Outlog.write("\t Number of Frameshift InDels selected (Med/High/Total): "+str(FrameshiftCount[0])+"/"+str(FrameshiftCount[1])+"/"+str(sum(FrameshiftCount))+"\n")
Outlog.write("\t Number of Splice site variants selected (Med/High/Total): "+str(SplicingCount[0])+"/"+str(SplicingCount[1])+"/"+str(sum(SplicingCount))+"\n")
Outlog.write("\t Number of Unknown function variants selected (Med/High/Total): "+str(UnknownCount[0])+"/"+str(UnknownCount[1])+"/"+str(sum(UnknownCount))+"\n")
Outlog.write("\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Filtering Finished: "+TimeNow+"\n")
print "Filtering finished "+TimeNow
