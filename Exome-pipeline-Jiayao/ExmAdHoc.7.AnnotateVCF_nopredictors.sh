#!/bin/bash
#$ -cwd -l mem=8G,time=12:: -N AnnVCF

#Same as pipeline script but does not added functional predictors such as Cadd, PP2 or snpEFF
#    InpFil - (required) - Path to VCF file to be annotated. Alternatively a file with a list of vcfs can be provided and the script run as an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $ANNHDB - directory containing databases for annovar

#list of required tools:
# annovar <http://www.openbioinformatics.org/annovar/> <http://www.openbioinformatics.org/annovar/annovar_download_form.php>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmVC.5.AnnotatewithANNOVAR.sh -i <InputFile> -r <reference_file> -l <logfile> -PH

     -i (required) - Path to VCF file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -H (flag) - echo this message and exit
"

PipeLine="false"
FullCadd="false"

while getopts i:r:l:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
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

#Set local Variables
##Set local parameters
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
AnnNam=`basename $InpFil`
AnnNam=${AnnNam/.vcf/}
if [[ -z $LogFil ]]; then LogFil=$AnnNam.AnnVCF.log; fi # a name for the log file
TmpLog=$AnnNam.AnnVCF.temp.log #temporary log file
AnnDir=$AnnNam.AnnVCF.tempdir; mkdir -p $AnnDir
TmpVar=$AnnDir/$AnnNam.tempvar
AnnFil=$AnnNam.annovar
VcfFil=${InpFil/vcf/annotated.vcf} # final annotated output file
TmpDir=$AnnNam.AnnVCF.tempdir; mkdir -p $TmpDir #temporary directory


#Start Log File
ProcessName="Annotate VCF" # Description of the script - used in log
funcWriteStartLog

##Convert VCF to ANNOVAR input file using ANNOVAR - use old vcf method
StepNam="Convert VCF to ANNOVAR input file using ANNOVAR"
StepCmd="convert2annovar.pl $InpFil -format vcf4old -includeinfo | cut -f 1-10 > $TmpVar" #command to be run
funcRunStep

##Generate Annotation table
StepNam="Build Annotation table using ANNOVAR"
StepCmd="table_annovar.pl $TmpVar $ANNHDB --buildver hg19 --remove -protocol refGene,esp6500si_all,esp6500si_aa,esp6500si_ea,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr, -operation g,f,f,f,f,f,f,f,f -otherinfo  -nastring \"\"  --outfile $AnnFil"
funcRunStep
AnnFil=$AnnFil.hg19_multianno.txt

##sort, replace spaces and semi-colons, zip and index
# - annovar output has spaces in the RefSeq function code, e.g. "synonymous SNV", but they are not permitted in vcf format and other tools (e.g. GATK) will throw an error if they encounter them
# - annovar separates multiple gene names in the RefSeq gene name field with semi-colons, this causes and error in the vcf
head -n 1 $AnnFil > $AnnFil.tempheader
tail -n+2 $AnnFil | sort -V | awk '{gsub( / /, ""); print}' | awk '{gsub( /;/, ","); print}' >> $AnnFil.tempheader
mv $AnnFil.tempheader $AnnFil
bgzip $AnnFil
tabix -S 1 -s 1 -b 2 -e 3 $AnnFil.gz

#Incorporate annovar annotations into vcf with vcftools
StepNam="Incorporate annovar annotations into vcf with vcftools"
StepCmd="cat $InpFil | vcf-annotate -a $AnnFil.gz
 -c -,-,-,-,-,INFO/VarFunc,INFO/GeneName,INFO/VarClass,INFO/AAChange,INFO/ESPfreq,-,-,INFO/1KGfreq,-,-,-,-,CHROM,POS,-,REF,ALT
 -d key=INFO,ID=VarFunc,Number=1,Type=String,Description='Genomic region/Sequence Function'
 -d key=INFO,ID=GeneName,Number=1,Type=String,Description='refGene GeneName'
 -d key=INFO,ID=VarClass,Number=1,Type=String,Description='Mutational Class'
 -d key=INFO,ID=AAChange,Number=1,Type=String,Description='Amino Acid change'
 -d key=INFO,ID=ESPfreq,Number=1,Type=Float,Description='Exome Sequencing Project 6500 alternative allele frequency'
 -d key=INFO,ID=1KGfreq,Number=1,Type=Float,Description='1000 genome alternative allele frequency'
 > $VcfFil"
funcRunStep

#End Log
funcWriteEndLog

#clean up
rm -rf $AnnDir
