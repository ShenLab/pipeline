
PipePath="$HOME/CUMC/Exome-Annotation/"
source $PipePath/mappability_generate.sh
function annotatedwithMappability {
   len=$1
   fMAP=$2
   fvcf=$3
   fbname=$(basename $fvcf)
   fbname=${fbname%.vcf.gz}
   fbname=${fbname%.vcf}
   fout="${fbname}.mappability.vcf.gz"
   echo $fMAP

  if [ ! -f $fMAP ]; then
         Redo_mappability $len
   fi
   keywd="Mappability"
   des="Mappability is obtained in file $fMAP"
   fhead="header_add.txt"

   echo "##INFO=<ID=$keywd,Number=1,Type=String,Description=\"$des\">" > $fhead
   bcftools annotate -a $fMAP  -c CHROM,FROM,TO,Mappability  -h  $fhead  $fvcf |bgzip -c > $fout
}


echo "input format annotatedwithMappability len *bed.gz  *.vcf.gz"
echo " bash $PipePath/mappability.annotation.sh len /path/to/hg19_lenbp_mappability.bed.gz  VCF.vcf"
annotatedwithMappability   $1 $2 $3
if [ ! -f $2 ];then
   echo "$2 cannot find!\n"
fi

if [ ! -f $3 ]; then
  echo "$3 cannot find!\n"
fi
