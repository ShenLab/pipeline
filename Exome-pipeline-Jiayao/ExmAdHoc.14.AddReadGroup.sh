#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N AddReadGroup 
#$ -l h_rt=24:00:00
#$ -l h_vmem=16G
#$ -cwd

# This script add or change the Read Groups for bam files

#set default arguments
usage="
-t 1-{number of fastq files] ExmAln.1a.Align_Fastq_to_Bam_with_BWAmem.sh -i <InputFile> -r <reference_file> -t <target intervals file> -l <logfile> -PH

     -i (required) - Table containing the path to the fastq file and the RG read header
     -r (required) - shell file containing variables with locations of reference files and resource directories (WES_Pipeline_References.b37.sh)
     -l (optional) - Log file
     -t (optional) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability); this file is required if calling the pipeline but otherwise can be omitted
     -P (flag) - Initiate exome analysis pipeline after completion of script
     -F (flag) - Fix mis-encoded base quality scores - Rewrite the qulaity scores from Illumina 1.5+ to Illumina 1.8+
     -H (flag) - echo this message and exit
"

PipeLine="false"

#get arguments
while getopts i:r:a:g:l:t:PFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        a) ArrNum="$OPTARG";; 
	    g) ReadGroup="$OPTARG";;
		l) LogFil="$OPTARG";;
        t) TgtBed="$OPTARG";; 
        P) PipeLine="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
    esac
done
echo "ReadGroup $ReadGroup "
echo "InpFil $InpFil "
echo "RefFil $RefFil "
#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi



#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil 

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#set local variables
#ArrNum=1 #line of table to read

InpFil=`readlink -f $InpFil`  # resolve input file path
BAM=`readlink -f $(tail -n+$ArrNum $InpFil | head -n 1 | cut -f1)`
BamNam=$(basename $BAM | sed s/.bam// ) # a name for the output files - basically the original file name


if [[ ! -e "$ReadGroup" ]];then
	ReadGroup=$(basename $BamNam|cut -d '.' -f 1)
fi
echo $ReadGroup
if [[ -z "$LogFil" ]]; then LogFil=$BamNam.FqB.log; fi # a name for the log file
ADRFil=$BamNam.tmp.bam
SrtFil=$BamNam.AddRG.bam #output file for sorted bam
DdpFil=$BamNam.mkdup.bam #output file for dedup bam
TmpLog=$BamNam.FqB.temp.log #temporary log file
TmpDir=$BamNam.FqB.tempdir; mkdir -p $TmpDir #temporary directory

#start log
ProcessName="Add Read Groups with Picard"
funcWriteStartLog
echo " Build of reference files: "$BUILD >> $TmpLog
echo "----------------------------------------------------------------" >> $TmpLog

StepName="Add Read Groups with Picard"
StepCmd="java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $PICARD AddOrReplaceReadGroups
 INPUT=$BAM
 OUTPUT=$ADRFil
 RGID=$ReadGroup
 RGLB=$ReadGroup
 RGPL=$ReadGroup
 RGPU=$ReadGroup
 RGSM=$ReadGroup
 2>>$TmpLog"
funcRunStep

#Sort the bam file by coordinate
StepName="Sort Bam using PICARD"
StepCmd="java -Xmx8G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=$TmpDir -jar $PICARD SortSam
 INPUT=$ADRFil
 OUTPUT=$SrtFil
 SORT_ORDER=coordinate
 CREATE_INDEX=TRUE
 2>>$TmpLog"
funcRunStep
rm $ADRFil #remove the "Aligned bam"

#Mark the duplicates
StepName="Mark PCR Duplicates using PICARD"
StepCmd="java -Xmx8G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=$TmpDir -jar $PICARD MarkDuplicates
 INPUT=$SrtFil
 OUTPUT=$DdpFil
 METRICS_FILE=$DdpFil.dup.metrics.txt
 CREATE_INDEX=TRUE
 2>>$TmpLog"
funcRunStep
rm $SrtFil ${SrtFil/bam/bai} #remove the "Sorted bam"

#End Log
funcWriteEndLog
echo $DdpFil
