#Pipeline for single family with PSAP
while getopts i:v:a:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        v) VCFFil="$OPTARG";;
        a) ArrNum="$OPTARG";;
        H) echo "$usage"; exit;;
  esac 
done


if [[ ! -e "$InpFil" ]] || [[ ! -e "$VCFFil" ]] ; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

InpFil=`readlink -f $InpFil` #resolve absolute path to bam
Indv=$(tail -n+$ArrNum $InpFil | head -n 1)
OutName=`basename $Indv | sed s/.list//`
#WorkDir=`dirname $Indv`
LogFil=$OutName.log
#cd $WorkDir

step1="vcftools --gzvcf $VCFFil --recode --recode-INFO-all --out $OutName --keep $Indv"
echo $step1 >$LogFil
$step1 >>$LogFil 2>>$LogFil

step2="python $HOME/CUMC/Exome-Filters-Jiayao/Clean.py -v $OutName.recode.vcf -o $OutName.vcf" 
echo $step2 >>$LogFil
$step2 >>$LogFil 2>>$LogFil
rm $OutName.recode.vcf

step3="python $HOME/CUMC/Exome-Filters-Jiayao/ExmFilter.py -v $OutName.vcf -p $OutName.ped -o $OutName.RareCoding.vcf"
$step3 >>$LogFil 2>>$LogFil
Signal=`tail -n 1 $LogFil`
if [ $Signal != "Done" ];
then
	echo "Step3 not successful." >> $LogFil
	exit
fi
rm $OutName.vcf

if [ -f $OutName.avinput ];
then
	rm -rf annotated $Outname.avinput
fi
step4="$HOME/CUMC/psap/family_psap_pipeline.sh $OutName.RareCoding.vcf $OutName $OutName.ped"
echo $step4 >>$LogFil
$step4 >>$LogFil 2>>$LogFil


step5="PSAP_Report.py -r ./annotated/$OutName.report.txt -p $OutName.ped -v $OutName.RareCoding.vcf"
echo $step5 >>$LogFil
$step5 >>$LogFil 2>>$LogFil

