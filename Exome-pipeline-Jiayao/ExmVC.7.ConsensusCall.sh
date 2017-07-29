#!/bin/bash

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $DBSNP - dbSNP vcf from GATK
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
-t 1-<X> ExmVC.2.GenotypeGVCFs.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

-i (required) - List of BAM files. List file name must end \".list\"
-r (required) - shell file containing variables with locations of reference files and resource directories
-t (required) - Exome capture kit targets or other genomic intervals bed file. 
-n (optional) - Analysis/output VCF name - will be derived from input filename if not provided; only used if calling pipeline
-l (optional) - Log file
-j (optional) - number of jobs in array; normally this is derived automatically from the -t argument of qsub command; this can be used to rerun individual sections that failed
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -X (flag) - Do not run Variant Quality Score Recalibration - only if calling pipeline
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

while getopts i:r:a:t:n:l:j:PXBH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";;
        a) ArrNum="$OPTARG";; 
        t) TgtBed="$OPTARG";;
        n) VcfNam="$OPTARG";;
        j) JobNum="$OPTARG";;
        l) LogFil="$OPTARG";;
        P) PipeLine="true";;
        X) NoRecal="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]] || [[ -z "$TgtBed" ]]; then
 echo "Missing/Incorrect required arguments"
 echo "provided arguments: -i $InpFil -r $RefFil -t $TgtBed"
 echo "usage: $usage"
 exit
fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"


#Set local Variables
funcGetTargetFile #If the target file has been specified using a code, get the full path from the exported variable
# The target file needs to be divided evenly between all the jobs. i.e. if the target file is 1000 lines long and there are 40 jobs, each job should have 25 lines of the target file
# bash arithmetic division actually gives the quotient, so if there are 1010 lines and 40 jobs the division would still give 25 lines per a job and the last 10 lines would be lost
# to compensate for this we will find the remainder (RemTar) and then add an extra line to the first $RemTar jobs

if [[ $JobNum ]]; then NumJobs=$JobNum; fi
echo $NumJobs
TarLen=$(cat $TgtBed | wc -l) 
RemTar=$(( TarLen % NumJobs )) # get remainder of target file length and number of jobs
QuoTar=$(( TarLen / NumJobs )) # get quotient of target file length and number of jobs
SttLn=1
DivLen=0
echo $RemTar
echo $SttLn
for ((i=1; i <= $ArrNum; i++)); do
    SttLn=$(( SttLn + DivLen ))
    if [[ $i -le $RemTar ]]; then
        DivLen=$(( QuoTar + 1 ))
        else
        DivLen=$QuoTar
    fi
done


if [[ -z "$VcfNam" ]];then VcfNam=`basename $InpFil`; VcfNam=${VcfNam/.list/}; fi # a name for the output files
if [[ -z $LogFil ]]; then LogFil=$VcfNam.GgVCF.log; fi # a name for the log file
VcfDir=$VcfNam.splitfiles; mkdir -p $VcfDir # Directory to output slices to
PrgDir=$VcfNam.progfiles; mkdir -p $PrgDir # "Progress directory" to output completion logs to, these are used by the merge script to check that all jobs in this array have finished (if jobs run out of time the hold is released even if the jobs did not complete)
VcfNam=$VcfNam.$ArrNum
PrgFil=$VcfNam.genotypingcomplete
VcfFil=$VcfDir/$VcfNam.vcf #Output File
VcfAnnFil=$VcfDir/$VcfNam.ann.vcf
VcfLeftAlnFil=$VcfDir/$VcfNam.LA.vcf
GatkLog=$VcfNam.GgVCF.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$VcfNam.GgVCF.temp.log #temporary log file
TmpDir=$VcfNam.GgVCF.tempdir; mkdir -p $TmpDir #temporary directory
TgtFil=$TmpDir/Range.$VcfNam.bed #exome capture range
tail -n+$SttLn $TgtBed | head -n $DivLen > $TgtFil #get exome capture range
infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A FisherStrand -A InbreedingCoeff -A QualByDepth -A ChromosomeCounts -A GenotypeSummaries -A StrandOddsRatio -A DepthPerSampleHC"
# -A HomopolymerRun
# -A SpanningDeletions 

#Start Log File
ProcessName="Running Platypus Calling Variants on Give Bam List" # Description of the script - used in log
funcWriteStartLog
echo "Target file line range: $SttLn - $(( $SttLn + $DivLen - 1 ))" >> $TmpLog

##Run Joint Variant Calling
StepName="Joint call Variants with Platypus"
StepCmd="python $PLATYPUS
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

##Left Align variants
StepName="Left align variants in the VCF with GATK"
StepCmd="java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T LeftAlignAndTrimVariants
 -R $REF
 -V $VcfFil
 -o $VcfLeftAlnFil
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep
mv -f $VcfLeftAlnFil $VcfFil 
mv -f $VcfLeftAlnFil.idx $VcfFil.idx

##Write completion log
touch $PrgDir/$PrgFil


#End Log
funcWriteEndLog
