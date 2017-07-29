#!/bin/bash
#$ -l mem=4G,time=6:: -cwd -S /bin/bash -N TerModBam

# This script is to fix a problem with paired-end Bam files that were aligned with old verions of bwa in which "/3" or "/4" was added to the mate in each read pair rather than /1 and /2, as in:
#D8GSQ5P1:4:1102:10621:68803#0    73    1    10009    0    101M    =    10009    0    ACCCTAACCCTAACCCTAACCCTA....
#D8GSQ5P1:4:1102:10621:68803#0/3    133    1    10009    0    *    =    10009    0    GTTAGGGTTAGGGTTAGGGCTGGG....
# This causes problems with the bam2fq-->bwa pipeline - the bam2fq does not recoginise the pairs and doesn't interleave them properly and hence bwa then sees them as all singletons.
#    InpBam - (required) - The bam file to be modified
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    OutBam - (optional) - A name for the output file. If this is not provided provided \"fixed\" will be added into the original filename.
#    Help - H - (flag) - get usage information

#list of required tools:
# samtools <http://samtools.sourceforge.net/> <http://sourceforge.net/projects/samtools/files/>


## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAdHoc.4.TheTeresaMod.sh -i <InputName> -r <reference_file> -o <OutputName>

     -i (required) - Bam file to be modified
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -o (optional) - Output filename - if not provided \"fixed\" will be added into the original filename
     -H (flag) - echo this message and exit
"

#get arguments
while getopts i:o:r:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";;
        o) OutFil="$OPTARG";;
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

#check all required paramaters present
if [[ ! -e "$InpFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Set local Variables
LogFil=${InpFil/.bam/.TeresaMod.log}
TmpLog=$LogFil.TeresaMod.temp.log #temporary log file
if [[ -z $OutFil ]]; then OutFil=`basename $InpFil | sed 's/bam$/fixed.bam/'`; fi



#Start Log File
ProcessName="Modify Bam file to remove unwanted trailing characters" # Description of the script - used in log
funcWriteStartLog

##Run genomic VCF generation
StepNam="Modify bam file using samtools"
StepCmd="samtools view -h $InpFil | 
 awk 'BEGIN { OFS = \"\t\" } { gsub(/#0\/3\$/,\"#0\",\$1); gsub(/#0\/4\$/,\"#0\",\$1); print \$0 }' | 
 samtools view -bS - > $OutFil"
funcRunStep

#End Log
funcWriteEndLog

