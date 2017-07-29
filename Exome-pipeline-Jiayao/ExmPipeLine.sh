#!/bin/bash
RefFil=$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
FastQC=
BWA=
MergeBAM=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
DoC=
HapCaller=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAln.2.HaplotypeCaller_GVCFmode.sh
JointGT=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.1hc.GenotypeGVCFs.sh
VCFMerge=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
VQSR=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.4.RecalibrateVariantQuality.sh
Annovar=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/Test.AnnotateVCF_direct.sh

ProjectHome=/home/local/users/jw/Analysis/0613/
ProjectName=PAH_16_393914
BedFil=/home/local/users/jw/resources/references/b37/CaptureKitBeds/CCDS_hg19.bed
FastQList=
RawBamList=
BamList=${ProjectHome}/PAH.16_393914.bam.list 
GVCFList=${ProjectHome}PAH.16_393914.g.vcf.list
SplitedDir=${ProjectHome}/JointGenotyping/${ProjectName}.splitfiles
RawVCF=${ProjectHome}/JointGenotyping/${ProjectName}.rawvariants.vcf.gz
VQSRVCF=${ProjectHome}/JointGenotyping/${ProjectName}.rawvariants.recalibrated.vcf
AnnovarVCF=${ProjectHome}/JointGenotyping/${ProjectName}.rawvariants.recalibrated.vcf.hg19_multianno.vcf.gz

mkdir -p ${ProjectHome}/JointGenotyping

#HaploytypeCaller
#NJob=`wc -l $BamList|cut -f 1 -d ' '`
#echo $NJob
#seq $NJob | parallel -j 20 --eta $HapCaller -i $BamList -r $RefFil -t $BedFil -a {}

#Joint Genotyping
GVCFList=/home/local/users/jw/Analysis/0613/PAH.16_393914.g.vcf.list
cd ${ProjectHome}/JointGenotyping
#seq 20 | parallel -j 20 --eta sh $JointGT -i $GVCFList -r $RefFil -a {} -j 20 -t $BedFil -n $ProjectName

#MergeVCF
#echo $SplitedDir
#$VCFMerge -i $SplitedDir -r $RefFil 

#VQSR
#echo $RawVcf
#$VQSR -i $RawVCF -r $RefFil 

#Annovar
if [ -e $VQSRVCF ]
then
	$Annovar -i $VQSRVCF -r $RefFil 
else
	$Annovar -i $RawVCF -r $RefFil 
fi
