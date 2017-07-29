#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N MergeGVCF
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -cwd

#This script takes a list of gVCF files (filename must end ".list") and combines them into a single file
#    InpFil - (required) - A list of gVCFs to be combined. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    VcfNam - (optional) - A name for the analysis - to be used for naming output files. Will be derived from input filename if not provided; only used if calling pipeline
#    LogFil - (optional) - File for logging progress
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format
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
ExmAdHoc.6.Merge_gVCFs.sh -i <InputFile> -r <reference_file> -o <OutputName> -l <logfile> -PABH


     -i (required) - \".list\" file containing a paths to gVCFs to be merged
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -o (optional) - Analysis/output VCF name - will be derived from input filename if not provided
     -l (optional) - Log file
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

BadET="false"

while getopts i:r:o:l:PBH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        o) VcfNam="$OPTARG";;
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
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"


#Set local Variables
funcGetTargetFile #If the target file has been specified using a code, get the full path from the exported variable
if [[ -z "$VcfNam" ]];then VcfNam=`basename $InpFil`; VcfNam=${VcfNam%%.*}; fi  # a name for the output files
if [[ -z $LogFil ]]; then LogFil=$VcfNam.CmbGVCF.log; fi # a name for the log file
InpFil=`readlink -f $InpFil` #resolve absolute path to bam
VcfFil=$VcfNam.Combined.g.vcf #Output File
GatkLog=$VcfNam.CmbGVCF.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$VcfNam.CmbGVCF.temp.log #temporary log file
TmpDir=$VcfNam.CmbGVCF.tempdir; mkdir -p $TmpDir #temporary directory

#Start Log File
ProcessName="Merge Genomic VCF" # Description of the script - used in log
funcWriteStartLog

##Run genomic VCF generation
StepNam="gVCF generation with GATK HaplotypeCaller"
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T CombineGVCFs
 -R $REF
 -V $InpFil
 -o $VcfFil
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

##gzip and index the gVCF
StepName="gzip and index the gVCF"
StepCmd="bgzip $VcfFil; tabix -f -p vcf $VcfFil.gz"
funcRunStep
rm $VcfFil.idx

#End Log
funcWriteEndLog
