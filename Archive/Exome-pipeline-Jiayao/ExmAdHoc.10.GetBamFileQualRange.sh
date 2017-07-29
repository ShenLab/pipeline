#!/bin/bash
#$ -cwd -l mem=4G,time=:20: -N GetQual

# This script output the maximum and minimum base quality scores observed in a bam file. If the bamfile contains an OQ fields (old quality scores retained after BQSR with GATK) the script will output the results for that field
# Usage notes:
# Please see ExmAln.1b.ReAlign_Bam_with_BWAmem.sh for full details 
# The script should be run from within the mapping directory
#    InpFil - (required) - Path to the Bam file
#    RedNum - (optional) - Number of reads of Bam file to be checked (default=1000)
#    OutFil - (optional) - Output file name. Will output to std out if not given
#    Help - H - (flag) - get usage information

#list of required tools:
# samtools <http://samtools.sourceforge.net/> <http://sourceforge.net/projects/samtools/files/>

###############################################################

#set default arguments
usage="
ExmAdHoc.10.GetInpFileQualRange.sh -i <InputFile> -n <NumberofReads> -o <OutputFilename> -H

     -i (required) - Aligned bam file
     -n (required) - Number of reads of Bam file to be checked (default=1000)
     -o (optional) - Output file name
     -H (flag) - echo this message and exit
"

RedNum=1000
#get arguments
while getopts i:n:o:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";; 
        n) RedNum"$OPTARG";; 
        o) OutFil="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

QualField=11
SamRead=`samtools view $InpFil | head -n 1`
#check for OQ and change quality field number if necessary

if [[ "$SamRead" == *OQ:* ]]; then
    QualAscii=`samtools view $InpFil | head -n $RedNum | awk '{ gsub( /.*[[:space:]]OQ:.:/, ""); gsub ( /[[:space:]].*/, ""); print }' | grep -o . | LC_ALL=C sort | uniq`
else
    QualAscii=`samtools view $InpFil | head -n $RedNum | cut -f $QualField | grep -o . | LC_ALL=C sort | uniq`
fi
minQual=${QualAscii:0:1}
maxQual=${QualAscii: -1}

if [[ $OutFil ]]; then
    echo $InpFil": Min="$minQual"; Max="$maxQual >> $OutFil
fi
echo $InpFil": Min="$minQual"; Max="$maxQual
