#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N FilterPC

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts v:t:n: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        t) PartFil="$OPTARG";;
        n) DirPre="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`
PartFil=`readlink -f $PartFil`

FamNam=`cut -f 1 $PartFil | head -n $SGE_TASK_ID | tail -n 1`
Proband=`cut -f 2 $PartFil | head -n $SGE_TASK_ID | tail -n 1`
Parent=`cut -f 3 $PartFil | head -n $SGE_TASK_ID | tail -n 1`


echo $FamNam
if [[ $FamNam == [0-9]* ]]; then FamNam=Fam$FamNam; fi

DirNam=$FamNam.Trio
if [[ -n $DirPre ]]; then DirNam=$DirPre"_"$DirNam; fi
mkdir -p $DirNam
cd $DirNam

#Autosomal Recessive
echo "Autosomal Recessive.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.AR --alt $Proband --het $Parent"
echo $CMD
eval $CMD
#Autosomal Dominant
echo "Autosomal Dominant.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.AD  --het $Proband --ref $Parent"
echo $CMD
eval $CMD
#compound heterozygous
echo "Compund heterozygous.."
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.tempheppat  --het $Proband,$Parent"
echo $CMD
eval $CMD
CMD="$FiltScrDir/ExmFilt.CustomGenotype.py -v $VcfFil -o $FamNam.temphepmat  --het $Proband --ref $Parent"
echo $CMD
eval $CMD
R --vanilla <<RSCRIPT
options(stringsAsFactors=F)

mathet <- read.delim("$FamNam.temphepmat.tsv")
pathet <- read.delim("$FamNam.tempheppat.tsv")

mathet <- mathet[mathet[,"Gene"]%in%pathet[,"Gene"],]
pathet <- pathet[pathet[,"Gene"]%in%mathet[,"Gene"],]

comphet <- rbind(mathet, pathet)
comphet <- comphet[order(comphet[,"Chromosome"], comphet[,"Position"]),]
write.table(comphet, "$FamNam.compound_heterozygous.tsv", col.names=T, row.names=F, quote=F, sep="\t")
RSCRIPT
cat $FamNam.tempheppat.log $FamNam.temphepmat.log > $FamNam.filter.compound_heterozygous.log

rm -rf *temp*
