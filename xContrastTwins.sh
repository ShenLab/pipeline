#!/bin/bash
#$ -cwd -l mem=1G,time=:30: -N FilterTrio

#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi

#get arguments
while getopts v:a:b:u:o:: opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        a) Twin1="$OPTARG";;
        b) Twin2="$OPTARG";;
        u) Unfiltered="$OPTARG";;
        o) OutFil="$OPTARG";;
        #H) echo "$usage"; exit;;
    esac
done

FiltScrDir="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/"

VcfFil=`readlink -f $VcfFil`

PyCall="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts/ExmFilt.CustomGenotype.py -v $VcfFil -o $OutFil.temp -f 0.01 -p " 

runCustomFilter(){
if [[ $Unfiltered ]]; then TwinComp=$TwinComp" --unfl $Unfiltered"; fi
CMD="$PyCall $TwinComp"
echo $CMD
eval $CMD
nlines=$(cat $OutFil.temp.tsv | wc -l)
nlines=$(( $nlines - 1 ))
echo "Filtered $nlines variants"
}

TwinComp=" --het $Twin1 --ref $Twin2"
runCustomFilter
cat $OutFil.temp.tsv >> $OutFil.tsv

TwinComp=" --het $Twin1 --alt $Twin2"
runCustomFilter
cat $OutFil.temp.tsv | awk '{ out=$1; for(i=2;i<=15;i++){out=out"\t"$i}; out=out"\t"$17"\t"$16"\t"$18; print out}' >> $OutFil.tsv

TwinComp=" --ref $Twin1 --het $Twin2"
runCustomFilter
cat $OutFil.temp.tsv | awk '{ out=$1; for(i=2;i<=15;i++){out=out"\t"$i}; out=out"\t"$17"\t"$16"\t"$18; print out}' >> $OutFil.tsv

TwinComp=" --ref $Twin1 --alt $Twin2"
runCustomFilter
cat $OutFil.temp.tsv | awk '{ out=$1; for(i=2;i<=15;i++){out=out"\t"$i}; out=out"\t"$17"\t"$16"\t"$18; print out}'  >> $OutFil.tsv

TwinComp=" --alt $Twin1 --ref $Twin2"
runCustomFilter
cat $OutFil.temp.tsv >> $OutFil.tsv

TwinComp=" --alt $Twin1 --het $Twin2"
runCustomFilter
cat $OutFil.temp.tsv >> $OutFil.tsv
rm $OutFil.temp*
