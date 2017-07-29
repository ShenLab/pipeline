#!/bin/bash
#$ -cwd -l mem=10G,time=4:: -N VcfPCA

#This script takes VCF file and runs PCA on the variants using Eigenstrat.
#    InpFil - (required) - A vcf file or plink bed/bim/fam trio of files. If plink, supply the bed file name (e.g. myfile.bed)
#    OutNam - (optional) - A name for the analysis - to be used for naming output files. Will be derived from input filename if not provided
#    Help - H - (flag) - get usage information

#list of required tools:
# 


###############################################################

usage="
ExmAdHoc.5.VCF_PCA.sh -i <InputFile> -l <logfile> -H

     -i (required) - A vcf input file
     -o (optional) - output name - will be derived from input filename if not provided
     -H (flag) - echo this message and exit
"

SamOnly="false"
while getopts i:o:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        o) OutNam="$OPTARG";;
        H) echo "$usage"; exit;;
  esac
done

#some variables
EXOMFILT=~/scripts/Filtering_scripts
HapMapReference=~/resources/hapmap3_pop/hg19/All_HapMap #root nhame for bed/bim/fam of HapMap data


InpFil=`readlink -f $InpFil`
if [[ -z "$OutNam" ]];then OutNam=`basename $InpFil`; OutNam=${OutNam/.bed/}; OutNam=${OutNam/.vcf/}; fi # a name for the output files
LogFil=$OutNam.PCA.log


#check for vcf, if vcf convert to plink format
if [[ "${InpFil##*.}" != "bed" ]]; then
    VcfFil=`readlink -f $InpFil`
    BbfNam=$OutNam
    #filter the VCF for common variants for 1KG and GO-ESP and within cohort alternate allele frequency > 0.05
    #$EXOMFILT/ExmFilt.1.FilterbyAlleleFrequency.py -v $VcfFil -o $BbfNam --maf 0.05 -G -W
    VcfFil=$BbfNam.filter.aaf.vcf
    #bgzip $BbfNam.filter.aaf.vcf.gz
    VcfFil=$BbfNam.filter.aaf.vcf.gz
    tabix -p vcf $VcfFil
    # convert vcf --> plink using vcftools; if there are more than 1000 samples in the vcf vcftools cannot generate the necesary number of temporary files (cluster limitiation) so we need to do it in multiple batches and then remerge
    SamLst=$OutNam.samples.list
    less $VcfFil | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' > $SamLst
    SamLen=`cat $SamLst | wc -l`
    if [[ $SamLen -lt 500 ]]; then
        vcftools --gzvcf $VcfFil --plink --out $BbfNam.tmp
        plink --file $BbfNam.tmp --make-bed --out $BbfNam
        rm -f $BbfNam.tmp*
    else
        mkdir -p TempSplit
        split --additional-suffix=.list -l 500 $SamLst $OutNam.samples.split.
        #for i in $OutNam.samples.split.*.list; do
        #    nohup vcftools --gzvcf $VcfFil --keep $i --plink --out TempSplit/${i/.list/} &
        #done
        #wait
        ls TempSplit/*ped | awk -v OFS='\t' '{ $2=$1 ; gsub ( /ped$/, "map", $2) ; print }' > $BbfNam.tmp.splitlist
        plink --merge-list $BbfNam.tmp.splitlist --make-bed --out $BbfNam
        #rm -rf TempSplit
    fi
    echo
    echo "------------------------------------------------------------------------"
    echo
    #change -9 in the fam to 2
    awk '{ gsub( /-9$/, "2"); print }' $BbfNam.fam > $BbfNam.fam.temp
    mv -f $BbfNam.fam.temp $BbfNam.fam
else
    BbfNam=`readlink -f $InpFil`
    BbfNam=${BbfNam/.bed/}
    awk '{ gsub( /-9$/, "2"); print }' $BbfNam.fam > $BbfNam.fam.temp
    mv $BbfNam.fam $BbfNam.fam.pcabkp
    mv $BbfNam.fam.temp $BbfNam.fam
fi

#Get a list of SNPs to retrieve from the reference data
SnpList=$OutNam.snplist
cut -f 2 $BbfNam.bim > $SnpList

# Get HapMap data
HapMapDat=$OutNam"_HapMapData"
plink --bfile $HapMapReference --extract $SnpList --allow-no-sex --make-bed --out $HapMapDat
if [[ $? -ne 0 ]]; then exit; fi
echo
echo "------------------------------------------------------------------------"
echo

###HapMap data is b36 so update map before merging
cut -f 2,4 $BbfNam.bim > update_map.tab
plink --bfile $HapMapDat --update-map update_map.tab --allow-no-sex --make-bed --out $HapMapDat
if [[ $? -ne 0 ]]; then exit; fi
echo
echo "------------------------------------------------------------------------"
echo

#change -9 in fam to 1 
awk '{ gsub( /-9$/, "1"); print }' $HapMapDat.fam > $HapMapDat.fam.temp
mv -f $HapMapDat.fam.temp $HapMapDat.fam

# Merge HapMap data:
EigDat=$OutNam.HapMap
plink --bfile $BbfNam --bmerge $HapMapDat --geno 0.05 --allow-no-sex --make-bed --out $EigDat
if [[ ! -e $EigDat-merge.missnp && $? -ne 0 ]]; then exit; fi
echo
echo "------------------------------------------------------------------------"
echo

#check for mismatched snps and multiple position/chr snps and exclude and remerge if necessary
grep "Warning: Multiple [cp]" $EigDat.log | sed s/.*\ \'//g | sed s/\'.*//g > ExcludeSNPs.list
cat $EigDat-merge.missnp >> ExcludeSNPs.list

if [[ -s ExcludeSNPs.list ]]; then
    plink --bfile $HapMapDat --exclude ExcludeSNPs.list --allow-no-sex --make-bed --out $HapMapDat
    if [[ $? -ne 0 ]]; then exit; fi
    echo
    echo "------------------------------------------------------------------------"
    echo
    plink --bfile $BbfNam --bmerge $HapMapDat --geno 0.05 --allow-no-sex --make-bed --out $EigDat
    if [[ ! -e $EigDat-merge.missnp && $? -ne 0 ]]; then exit; fi
    echo
    echo "------------------------------------------------------------------------"
    echo
fi

# Convert data to Eigenstrat format
cp $EigDat.fam $EigDat.pedind
echo genotypename: $EigDat.bed > par.BBF.EIGENSTRAT
echo snpname: $EigDat.bim >> par.BBF.EIGENSTRAT
echo indivname: $EigDat.pedind >> par.BBF.EIGENSTRAT
echo outputformat: EIGENSTRAT >> par.BBF.EIGENSTRAT
echo genooutfilename: $OutNam.eigenstratgeno >> par.BBF.EIGENSTRAT
echo snpoutfilename: $OutNam.snp >> par.BBF.EIGENSTRAT
echo indoutfilename: $OutNam.ind >> par.BBF.EIGENSTRAT
echo "par.BBF.EIGENSTRAT:"
cat par.BBF.EIGENSTRAT
echo
echo "------------------------------------------------------------------------"
echo

convertf -p par.BBF.EIGENSTRAT

if [[ $? -ne 0 ]]; then exit; fi
echo
echo "------------------------------------------------------------------------"
echo

# run EigenStrat

CMD="smartpca.perl -i $OutNam.eigenstratgeno -a $OutNam.snp -b $OutNam.ind -k 10 -o $OutNam.plus.HapMap.pca -p $OutNam.plus.HapMap.plot -e $OutNam.plus.HapMap.eval -l $OutNam.plus.HapMap.log -m 5 -t 2 -s 6.0"
echo $MCD
eval $CMD

if [[ -e $BbfNam.fam.pcabkp ]]; then mv $BbfNam.fam.pcabkp $BbfNam.fam; fi
