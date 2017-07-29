#!/bin/bash
#$ -cwd -l mem=10G,time=2:: -N GATKAnn

#This script takes a vcf file and uses GATK to reannotate with various variant quality and calling information.
#    InpFil - (required) - A vcf file to be annotated
#    RefFil - (required) - shell file containing variables with locations of reference files and resource directories; see list below for required variables for required variables
#    LogFil - (optional) - File for logging progress
#    Help - H - (flag) - get usage information

#list of required variables in Reference File:
# $REF - reference genome in fasta format
# $DBSNP - dbSNP vcf from GATK
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAdHoc.3.GATKAnnotateVCF.sh -i <InputFile> -r <reference_file> -t <target intervals file> -l <logfile> -H

     -i (required) - A vcf file to be annotated
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

BadET="false"

#get arguments
while getopts i:r:l:BH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil 

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#set local variables
VcfFil=`readlink -f $InpFil` #resolve absolute path to Vcf
VcfNam=`basename $VcfFil | sed s/.vcf//` # a name for the output files
if [[ -z "$LogFil" ]]; then LogFil=$VcfNam.AnnGATK.log; fi # a name for the log file
VcfAnnFil=$VcfNam.annotated.vcf #output file with PCR duplicates marked
GatkLog=$VcfNam.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$VcfNam.AnnGATK.temp.log #temporary log file
TmpDir=$VcfNam.AnnGATK.tempdir; mkdir -p $TmpDir #temporary directory
infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A HaplotypeScore -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions -A FisherStrand -A InbreedingCoeff" #Annotation fields to output into vcf files

#start log
ProcessName="Merge with samtools"
funcWriteStartLog

##Annotate VCF with GATK
StepNam="Joint call gVCFs" >> $TmpLog
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T VariantAnnotator 
 -R $REF
 -L $VcfFil
 -V $VcfFil
 -o $VcfAnnFil
 -D $DBSNP
  $infofields
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#End Log
funcWriteEndLog
