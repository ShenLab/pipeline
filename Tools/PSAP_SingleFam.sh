#Pipeline for single family with PSAP
while getopts p:v:H opt; do
    case "$opt" in
        p) PedFil="$OPTARG";;
        v) VCFFil="$OPTARG";;
        H) echo "$usage"; exit;;
  esac 
done


if [[ ! -e "$VCFFil" ]] || [[ ! -e "$PedFil" ]] ; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi
PedFil=`readlink -f ${PedFil}`
FAMID=$(basename $PedFil| sed s/.ped//)
ListFil=${FAMID}.list
VCFFil=`readlink -f ${VCFFil}`
LogFil=$FAMID.log
cat $PedFil |grep -v "^#" |cut -f2 > $ListFil 
step1="vcftools --gzvcf $VCFFil --recode --recode-INFO-all --out $FAMID --keep $ListFil"
echo $step1 >$LogFil
$step1 >>$LogFil 2>>$LogFil

step2="python $HOME/CUMC/Exome-Filters-Jiayao/Clean.py -v $FAMID.recode.vcf -o $FAMID.vcf" 
echo $step2 >>$LogFil
$step2 >>$LogFil 2>>$LogFil
rm $FAMID.recode.vcf

step3="python $HOME/CUMC/Exome-Filters-Jiayao/ExmFilter.py -v $FAMID.vcf -p $FAMID.ped -o $FAMID.RareCoding.vcf -f ${HOME}/CUMC/Exome-Filters-Jiayao/ALL_FILTER.yml"
$step3 >>$LogFil 2>>$LogFil
Signal=`tail -n 1 $LogFil`
if [[ $Signal != "Done" ]] ;
then
	echo "Step3 not successful." >> $LogFil
	exit
fi
#rm $FAMID.vcf

if [ -f $FAMID.avinput ];
then
	rm -rf annotated $Outname.avinput
fi
step4="$HOME/CUMC/psap/family_psap_pipeline.sh $FAMID.RareCoding.vcf $FAMID $FAMID.ped"
echo $step4 >>$LogFil
$step4 >>$LogFil 2>>$LogFil


step5="python $HOME/CUMC/Tools/PSAP_Report.py -r ./annotated/$FAMID.report.txt -p $FAMID.ped -v $FAMID.RareCoding.vcf"
echo $step5 >>$LogFil
$step5 >>$LogFil 2>>$LogFil

