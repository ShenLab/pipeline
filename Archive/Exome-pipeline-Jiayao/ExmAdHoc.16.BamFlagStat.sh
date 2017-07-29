#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N BamFlagStat 
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -cwd


#This script takes a bam file and generates depth of coverage statistics using GATK
#    InpFil - (required) - Path to Bam file to be aligned or a file containing a list of bam files one per line (file names must end ".list")
#            if it is a list then call the job as an array job with -t 1:n where n is the number of bams
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
#    LogFil - (optional) - File for logging progress
#    Flag - B - BadET - prevent GATK from phoning home
#    Flag - F - Fix mis-encoded base quality scores - see GATK manual. GATK will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $TARGET - exome capture intervals bed file or other target file (must end ".bed")
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
 (-t <X>-<Y> [if providing a list]) ExmAln.8a.DepthofCoverage.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -DCBH

     -i (required) - Path to Bam file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -t (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
     -l (optional) - Log file
     -D (flag) - keep full Depth of Coverage file
     -C (flag) - Allow BadCigar - see GATK documentation - allows reads that GATK interprets as indicating a malformed file, e.g. reads starting with a deletion
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

BadCigar="false"
BadEt="false"
FullDoC="false"

#get arguments
while getopts i:r:a:t:l:DCBFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";;
		a) ArrNum="$OPTARG";; 
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] ; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

RefFil=`readlink -f $RefFil`
source $RefFil
source $EXOMPPLN/exome.lib.sh
funcGetTargetFile

#Set Local Variables
InpFil=`readlink -f $InpFil` #resolve absolute path to bam
BamFil=$(tail -n+$ArrNum $InpFil | head -n 1) 
BamNam=`basename $BamFil | sed s/.bam//` #a name to use for the various files

#Start Log
ProcessName="Depth of Coverage with GATK" # Description of the script - used in log
funcWriteStartLog

#Calculate depth of coverage statistics
StepName="Samtools flagstat" # Description of this step - used in log
StepCmd="samtools flagstat $BamFil  >> $BamNam.reads.mapped" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#End Log
funcWriteEndLog
