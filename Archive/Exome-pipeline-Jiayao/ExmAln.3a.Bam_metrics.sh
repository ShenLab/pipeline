#!/bin/bash
#$ -cwd -l mem=8G,time=6:: -N BamMtr 

# This script takes a bam file and generates insert size, GC content and quality score metrics using Picard
#    InpFil - i - (required) - Path to Bam file.
#    RefFiles - r - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - l - (optional) - File for logging progress
#    Metrics: G, I, Q - (flags) - will run GC bias (G), Insert Size (I) or Quality Distribution (Q); default is to run all metrics, specifying one or more will only run those specified
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $PICARD - directory containing Picard jar files

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# picard <http://picard.sourceforge.net/> <http://sourceforge.net/projects/picard/files/>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################


usage="
ExmAln.3a.Bam_metrics.sh -i <InputFile> -r <reference_file> -l <logfile> -GIQH

     -i (required) - Path to Bam file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -G (flag) - get GC bias metrics
     -I (flag) - get Inset Size metrics
     -Q (flag) - get quality distribution metrics  **All three metirics are run be default
     -H (flag) - echo this message and exit
"

#get arguments
while getopts i:r:l:GIQH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";;
        l) LogFil="$OPTARG";;
        G) GCmet="true";;
        I) ISmet="true";;
        Q) QDmet="true";;
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

#Set local Variables
if [[ -z $GCmet ]] && [[ -z $ISmet ]] && [[ -z $QDmet ]]; then #if no flags run all metrics
    ALLmet="true"
fi
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil | sed s/.bam$//` #a name to use for the various output files
if [[ -z $LogFil ]];then
    LogFil=$BamNam.BamMetrics.log # a name for the log file
fi
TmpLog=$BamNam.BamMet.temp.log #temporary log file 
TmpDir=$BamNam.BamMet.tempdir; mkdir -p $TmpDir #temporary directory

#Start Log
ProcessName="Start Get GC metrics with Picard" # Description of the script - used in log
funcWriteStartLog

#Get GC metrics with Picard
if [[ $ALLmet == "true" ]] || [[ $GCmet == "true" ]]; then
    StepName="Get GC Metrics with Picard" # Description of this step - used in log
    StepCmd="java -Xmx8G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $PICARD CollectGcBiasMetrics
 INPUT=$BamFil
 OUTPUT=$BamNam.GCbias_detail
 SUMMARY_OUTPUT=$BamNam.GCbias_summary
 CHART=$BamNam.GCbias.pdf
 REFERENCE_SEQUENCE=$REF
 VALIDATION_STRINGENCY=SILENT
 WINDOW_SIZE=200" #command to be run
    funcRunStep
fi

#Get Insert size metrics with Picard
if [[ $ALLmet == "true" ]] || [[ $ISmet == "true" ]]; then
    StepName="Get Insert Size Metrics with Picard" # Description of this step - used in log
    StepCmd="java -Xmx8G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $PICARD CollectInsertSizeMetrics
 INPUT=$BamFil
 OUTPUT=$BamNam.InsertSize_detail
 HISTOGRAM_FILE=$BamNam.InsertSize.pdf
 VALIDATION_STRINGENCY=SILENT" #command to be run
    funcRunStep
fi

#Quality Score Distribution
if [[ $ALLmet == "true" ]] || [[ $QDmet == "true" ]]; then
    StepName="Get Quality Score Distribution from BAM file using PICARD" # Description of this step - used in log
    StepCmd="java -Xmx8G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $PICARD QualityScoreDistribution
 INPUT=$BamFil
 OUTPUT=$BamNam.QualityDistr
 CHART_OUTPUT=$BamNam.QualityScoreDistr.pdf
 REFERENCE_SEQUENCE=$REF
 VALIDATION_STRINGENCY=SILENT" #command to be run
    funcRunStep
fi
#End Log
funcWriteEndLog
