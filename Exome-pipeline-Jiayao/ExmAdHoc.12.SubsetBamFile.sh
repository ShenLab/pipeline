#!/bin/bash
#$ -l mem=4G,time=6:: -cwd -S /bin/bash -N TrimBam

# This script is used to subset a bam file using locations in a bed file
#    InpBam - (required) - The bam file to be modified
#    TgtBed - (required) - input bed file with the required regions
#    OutBam - (optional) - A name for the output file. If this is not provided provided \"fixed\" will be added into the original filename.
#    Help - H - (flag) - get usage information

#list of required tools:
# samtools <http://samtools.sourceforge.net/> <http://sourceforge.net/projects/samtools/files/>


###############################################################

#set default arguments
usage="
ExmAdHoc.12.SubsetBamFile.sh -i <InputName> -b <TgtBede> -o <OutputName>

     -i (required) - Bam file to be modified
     -o (optional) - Output filename - if not provided \"fixed\" will be added into the original filename
     -H (flag) - echo this message and exit
"

#get arguments
while getopts i:t:o:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        t) TgtBed="$OPTARG";;
        o) OutBam="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpBam" ]] || [[ ! -e "$TgtBed" ]]; then 
    echo "Missing/Incorrect required arguments:"
    echo "   Input Bam: "$InpFil
    echo "   Input Bed: "$TgtBed
    echo "$usage"
    exit
fi

BamFil=`readlink -f $InpFil` #resolve absolute path to bam
TgtBed=`readlink -f $TgtBed` #resolve absolute path to bam

if [[ -z $OutBam ]]; then
    OutBam=`basename $BamFil | sed 's/bam$/trimmed.bam/'`
fi
OutBam=${OutBam%.bam}.bam
while [[ -e $OutBam ]]; do
    NUM=$(( NUM + 1 ))
    OutBam=${OutBam%.bam}
    OutBam=${OutBam%.[0-9]*}
    OutBam=$OutBam.$NUM.bam
done

LogFil=${OutBam/.bam/.log}
echo "Start Subset - $0:`date`" >> $LogFil
echo "Input Bam: "$BamFil >> $LogFil
echo "Target bed file: "$TgtBed >> $LogFil
echo "Output Bam: "$OutBam >> $LogFil

samtools view -b -L $TgtBed $BamFil > $OutBam.temp
mv $OutBam.temp $OutBam
samtools index $OutBam

echo "Completed $0:`date`" >> $LogFil
qstat -j $JOB_ID | grep -E "usage" >> $LogFil
