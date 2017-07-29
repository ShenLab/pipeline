#!/bin/bash
#$ -cwd -l mem=12G,time=4:: -N LocRln


#This script takes a bam file andperfomrs local indel realignment using GATK, the script can be run as an array across 24 Chromosomes by adding -t 1-24
#	InpFil - (required) - Path to Bam file to be realigned
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	LogFil - (optional) - File for logging progress
#	Flag - B - BadET - prevent GATK from phoning home

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $INDEL - Gold standard INDEL reference from GATK
# $INDEL1KG - INDEL reference from 1000 genomes
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java
# GATK

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
BadEt="false"

#get arguments
while getopts i:r:l:B opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
		l) LogFil="$OPTARG";;
		B) BadET="true";;
	esac
done

#load settings file
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set Local Variables
if [[ $SGE_TASK_ID != "undefined" ]]; then
    Chr=$SGE_TASK_ID
    Chr=${Chr/23/X}
    Chr=${Chr/24/Y}
    if [[ "$BUILD" = "hg19" ]]; then
        Chr=chr$Chr
    fi
    ChrNam=.$Chr
    RalDir=LocRealign.$BamNam/ #directory to collect individual chromosome realignments
    mkdir -p $RalDir
fi
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename ${BamFil/.bam/}` #a name to use for the various files
TmpDir=$BamNam$ChrNam.LocRealignjavdir #temp directory for java machine
mkdir -p $TmpDir
if [[ -z $LogFil ]];then
	LogFil=$BamNam.LocReal.log # a name for the log file
fi
TmpLog=$LogFil.LocReal$ChrNam.log #temporary log file 

RalLst=LocRealign.$BamNam.list #File listing paths to individual chromosome realignments
StatFil=RealignStatus$ChrNam.$JOB_ID.LocReal.stat #Status file to check if all chromosome are complete
TgtFil=$RalDir$BamNam$ChrNam.target_intervals.list #temporary file for target intervals for the CHR
realignedFile=$RalDir$BamNam.realigned$ChrNam.bam # the output - a realigned bam file for the CHR
GatkLog=$BamNam$ChrNam.LocRealign.gatklog #a log for GATK to output to, this is then trimmed and added to the script log

#Start Log
ProcessName="Start Local Realignment around InDels on Chromosome $Chr" # Description of the script - used in log
funcWriteStartLog

#Generate target file
StepName="Create target interval file using GATK RealignerTargetCreator" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T RealignerTargetCreator
 -R $REF
 -I $BamFil
 -known $INDEL
 -known $INDEL1KG
 -o $TgtFil
 -log $GatkLog"
if [[ $Chr ]]; then StepCmd=$StepCmd" -L $Chr"; fi
 #command to be run
funcGatkAddArguments
funcRunStep

#Realign InDels
StepName="Realign InDels file using GATK IndelRealigner" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T IndelRealigner
 -R $REF
 -I $BamFil
 -targetIntervals $TgtFil
 -known $INDEL
 -known $INDEL1KG
 -o $realignedFile
 -log $GatkLog"
if [[ $Chr ]]; then StepCmd=$StepCmd" -L $Chr"; fi
 #command to be run
funcGatkAddArguments
funcRunStep

#generate realigned file list
#find `pwd` | grep -E bam$ | grep $RalDir | sort -V > $RalLst

#End Log
funcWriteEndLog
#Clean up
#rm -r $TmpDir $TmpLog $TgtFil