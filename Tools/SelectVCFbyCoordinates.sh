usage="SelectVCFbyCoordinates.sh -i <InputFile> -o <OutputFile> -c <Chrom> -s <Start> -e <End>"

while getopts i:o:c:s:e:H opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		o) OutFil="$OPTARG";;
		c) Chrom="$OPTARG";;
		s) Start="$OPTARG";;
		e) End="$OPTARG";;
		H) echo "$usage"; exit;;
	esac
done
echo "-i $InpFil -o $OutFil -c $Chrom -s $Start -e $End	"

#if [[ ! -e "$InpFil" ]] || [[ ! -e "$OutFil" ]] || [[ ! -e "$Chrom" ]];  
#then 
#	echo "Missing/Incorrect required arguments"; 
#	echo "$usage"; 
#	exit; 
#fi

tabix -H $InpFil > $OutFil;
tabix $InpFil $Chrom:$Start-$End >> $OutFil;
