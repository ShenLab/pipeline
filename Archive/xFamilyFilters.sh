#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterFam

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts v:t:n: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) FamFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`
FamFil=`readlink -f $FamFil`

FamNam=`cut -f 1 $FamFil | head -n $SGE_TASK_ID | tail -n 1`
ModNam=`cut -f 2 $FamFil | head -n $SGE_TASK_ID | tail -n 1`
FilPrm=`cut -f 3 $FamFil | head -n $SGE_TASK_ID | tail -n 1`

if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam.Fam
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

if [[ $ModNam == *X* ]]; then FilPrm=$FilPrm" -X"; fi
if [[ $ModNam == *ovo* ]]; then FilPrm=$FilPrm" -D"; fi

echo "Filtering.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.$ModNam $FilPrm"
echo $CMD
eval $CMD

