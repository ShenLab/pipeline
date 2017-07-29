#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N Trimming 
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
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
        t) Threads="$OPTARG";; 
        P) PipeLine="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

if [[ ! -e "$Threads" ]]; then
	Threads=1
fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil 

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#set local variables
InpFil=`readlink -f $InpFil`  # resolve input file path
NCOL=$(head -n1 $InpFil | wc -w | cut -d" " -f1) #get number of columns in file to determine SE or PE
fastq1=`readlink -f $(tail -n+$ArrNum $InpFil | head -n 1 | cut -f1)` #R1 or SE fastq from first column
OutNam=$(echo $fastq1 | sed s/.fq.gz// | sed s/.fastq.gz//) # a name for the output files - basically the original file name
Out_Paired1=$OutNam.Paired1.fq.gz
Out_UnPaired1=$OutNam.UnPaired1.fq.gz
if [ $NCOL -eq 3 ]; then
    fastq2=`readlink -f $(tail -n+$ArrNum $InpFil | head -n 1 | cut -f3)` #R2 fastq from third column if present
	OutNam=$(echo $fastq1 | sed s/.fq.gz// | sed s/.fastq.gz//) # a name for the output files - basically the original file name
	Out_Paired2=$OutNam.Paired2.fq.gz
	Out_UnPaired2=$OutNam.UnPaired2.fq.gz
else
    fastq2=""
fi
rgheader=$(tail -n+$ArrNum $InpFil | head -n 1 | cut -f2) #RG header from second column

TmpDir=$OutNam.dir
TmpLog=$OutNam.log

if [ "$fastq2" != "" ];then
	StepName="Trimming reads with Adaptors And Low Base Q Paired End"
	echo $StepName
	StepCmd="java -Xmx4G -XX:ParallelGCThreads=$Threads -Djava.io.tmpdir=$TmpDir -jar $Trimmomatic 
		PE -threads $Threads
		-phred33
		$fastq1 $fastq2
		$Out_Paired1 $Out_UnPaired1
		$Out_Paired2 $Out_UnPaired2
		ILLUMINACLIP:$Adapter_TruSeq3_PE:2:30:10 
		LEADING:30 TRAILING:30 SLIDINGWINDOW:4:20 MINLEN:36
		2>>$TmpLog"
	funcRunStep
else
	StepName="Trimming reads with Adaptors And Low Base Q Single End"
	echo $StepName
	StepCmd="java -Xmx4G -XX:ParallelGCThreads=$Threads -Djava.io.tmpdir=$TmpDir -jar $Trimmomatic 
		SE -threads $Threads
		-phred33
		$fastq1 
		$Out_Paired1 
		ILLUMINACLIP:$Adapter_TruSeq3_SE:2:30:10 
		LEADING:30 TRAILING:30 SLIDINGWINDOW:4:20 MINLEN:36
		2>>$TmpLog"
	funcRunStep
fi
#End Log
funcWriteEndLog
echo "Done Trimming"
