#!/bin/bash
#$ -cwd -l mem=16G,time=6:: -N RenameSam

# This script can be used to rename a sample in a BAM file. It requires the a bam file and a new name for the sample. It will replace the sample name in all RG lines in the bam file.
# Usage notes:
# Please see ExmAln.1b.ReAlign_Bam_with_BWAmem.sh for full details 
# The script should be run from within the mapping directory
#    InpFil - (required) - Path to the Bam file
#    NewNam - (required) - New sample name to be inserted
#    OutFil - (optional) - Output file name
#    LogFil - (optional) - File for logging progress
#    OverWrite -O - (flag) - Overwrite the original bam file
#    Help - H - (flag) - get usage information

#list of required tools:
# samtools <http://samtools.sourceforge.net/> <http://sourceforge.net/projects/samtools/files/>

#list of required variables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAdHoc.9.Alignment_Finisher.sh -i <InputFile> -n <NewSampleNam> -o <OutputFilename> -l <logfile> -OH

     -i (required) - Aligned bam file
     -n (required) - New Sample Name
     -o (optional) - Output file name
     -O (optional) - Overwrite original file
     -l (optional) - Log file
     -H (flag) - echo this message and exit
"

OverWrite="false"
#get arguments
while getopts i:n:o:l:OH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";; 
        n) NewNam="$OPTARG";; 
        o) OutFil="true";;
        O) OverWrite="true";;
        l) LogFil="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ -z "$NewNam" ]]; then 
    echo "Missing/Incorrect required arguments"
    echo "$usage"
    echo "Provided arguments:"
    echo "    InpFil: "$InpFil"    NewNam: "$NewNam
    exit
fi

#set local variables
if [[ -z "$OutFil" ]]; then OutFil=`echo $InpFil | sed s/.bam$/.reheader.bam/`; fi # a name for the output file
if [[ -z "$LogFil" ]]; then LogFil=${OutFil%bam}log; fi # a name for the log file
echo "    InpFil: "$InpFil"    NewNam: "$NewNam"    OutFil: "$OutFil
TmpLog=$InpFil.rhcs.temp.log #temporary log file
TmpHead=$InpFil.rhcs.temp.header #temporary header file

#Export Header
StepName="Export header using Samtools"
StepCmd="samtools view -H $InpFil > $TmpHead"
funcRunStep
echo "Current read group header(s):" >> $TmpLog
grep ^@RG $TmpHead >> $TmpLog

#Replace Sample ID
StepName="Replace Sample ID using awk"
StepCmd="awk '{ if ( \$1 == \"@RG\" ) gsub ( /SM:[[:graph:]]*/ , \"SM:$NewNam\" ); print \$0 }' $TmpHead > $TmpHead.2;
 mv $TmpHead.2 $TmpHead"
funcRunStep
echo "New read group header(s):" >> $TmpLog
grep ^@RG $TmpHead >> $TmpLog

#Reheader the bam using samtools
StepName="Export header using Samtools"
StepCmd="samtools reheader $TmpHead $InpFil > $OutFil"
funcRunStep
echo "Check read group header(s) in new bam:" >> $TmpLog
samtools view -H $OutFil | grep ^@RG >> $TmpLog

if [[ "$OverWrite" == "true" ]]; then
    #replace original file
    StepName="Overwrite original file"
    StepCmd="mv -f $OutFil $InpFil"
    funcRunStep
fi

rm -f $TmpHead

#End Log
funcWriteEndLog
