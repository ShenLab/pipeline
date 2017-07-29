#!/bin/bash
#$ -cwd -l mem=4G,time=6:: -N AnnVCF

MPCFILE="/home/local/users/jw/resources/Annotations/MPC.hg19.Correct.txt.gz"

BCFTOOLS="/home/local/users/jw/software_packages/bcftools-1.4/bcftools"

#set default arguments
usage="
ExmVC.5.AnnotatewithANNOVAR.sh -i <InputFile> -r <reference_file> -l <logfile> -PH
     -i (required) - Path to VCF file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -C (flag) - Annotate with full CADD database
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -X (flag) - Do not run Variant Quality Score Recalibration
     -B (flag) - Prevent GATK from phoning home - only if calling pipeline
     -H (flag) - echo this message and exit
"

while getopts i:r:t:l:CXFPBH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        t) Threads="$OPTARG";;
        C) FullCadd="true";;
        P) PipeLine="true";;
        X) NoRecal="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

OutFil=$(basename $InpFil | sed s/.gz// | sed s/.vcf//)
OutFil=$OutFil.AddMPC.vcf.gz
#OutFil=$OutFil.AddMPC.vcf
HEADER=tmp.hdr
echo "##INFO=<ID=MPC,Number=1,Type=String,Description=\"MPC Score\">" > $HEADER
# Annotate from a tab-delimited file with six columns (the fifth is ignored),
# first indexing with tabix. The coordinates are 1-based.
#tabix -s1 -b2 -e2 annots.tab.gz
#bcftools annotate -a annots.tab.gz -h annots.hdr -c CHROM,POS,REF,ALT,-,TAG file.vcf

$BCFTOOLS annotate -a $MPCFILE -h $HEADER -c CHROM,POS,REF,ALT,MPC  $InpFil|bgzip -c > $OutFil
#tabix -p vcf $OutFil
#$BCFTOOLS annotate -a $MPCFILE -h $HEADER -c CHROM,POS,REF,ALT,MPC -O v $InpFil -o $OutFil
#StepCmd="$BCFTOOLS annotate -a $MPCFILE -h $HEADER -c CHROM,POS,REF,ALT,MPC $InpFil|bgzip -c $OutFil"
#$StepCmd
rm $HEADER
