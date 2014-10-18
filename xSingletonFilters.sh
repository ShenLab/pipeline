#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterSing

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts v:t:n: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) SingFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`
SingFil=`readlink -f $SingFil`

Proband=`cut -f 2 $SingFil | head -n $SGE_TASK_ID | tail -n 1`
FamNam=`cut -f 1 $SingFil | head -n $SGE_TASK_ID | tail -n 1`

if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam.Sgtn
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#Autosomal Recessive
echo "Autosomal Recessive.."
$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.AR --alt $Proband -P
#Autosomal Dominant
echo "Autosomal Dominant.."
$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.AD  --het $Proband -P
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

het <- read.delim("$FamNam.AD.tsv")

gens <- unique(het[duplicated(het[,"Gene"]),"Gene"])

comphet <- het[het[,"Gene"]%in%gens,]

write.table(comphet, "$FamNam.filter.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
