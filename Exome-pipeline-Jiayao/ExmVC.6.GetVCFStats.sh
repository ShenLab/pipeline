#!/bin/bash
#$ -cwd -l mem=4G,time=2:: -N VcfStats

#This script takes a raw VCF file and performs GATK's variant quality score recalibration
#    InpFil - (required) - Path to VCF file or a list of VCF Files to be recalibrated
#    LogFil - (optional) - File for logging progress
#    Help - H - (flag) - get usage information


#list of required tools:
# python

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="ExmVC.6.GetVCFStats.sh -i <InputFile> -l <logfile> -H

     -i (required) - Path to list of Bam files for variant calling
     -l (optional) - Log file
     -H (flag) - echo this message and exit
"

while getopts i:l:H opt; do
    AllOpts=$AllOpts" -$opt $OPTARG"
    case "$opt" in
        i) InpFil="$OPTARG";;
        l) LogFil="$OPTARG";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]]; then echo "Missing/Incorrect required arguments"; echo "Provided opts: "$AllOpts; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set local Variables
VcfFil=`readlink -f $InpFil` #resolve absolute path to vcf
VcfNam=`basename $VcfFil | sed s/.gz$// | sed s/.vcf$// ` #basename for outputs
StatFil=$VcfNam.stats.tsv

#Start Log File
ProcessName="Get VCF variant counts" # Description of the script - used in log
funcWriteStartLog

#Get VCF stats with python script
StepName="Get VCF stats"
StepCmd="python $EXOMPPLN/ExmPY.VCF_summary_Stats.py -v $VcfFil -o $StatFil"
funcRunStep

#End Log
funcWriteEndLog
