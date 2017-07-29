#!/bin/bash
# This script generates collapsed IGV plots given variant position and bam location 
# the plots will be in the same folder of the given variant file
# example: sh Generate_IGV_plots.sh -b example/bam.txt -v example/IGV_variants.txt
#

###############################################################
#set default arguments
usage="
     -v (required) - Table containing the path to the fastq file and the RG read header
     -b (required) - shell file containing variables with locations of reference files and resource directories (WES_Pipeline_References.b37.sh)
     -E (flag) - generate IGV plots in expand mode
     -H (flag) - echo this message and exit
"

#get arguments
ExpandMode="false"
while getopts v:b:g:EH opt; do
    case "$opt" in
        v) INDELS="$OPTARG";;
        b) BAM="$OPTARG";; 
		g) HG="$OPTARG";;
		E) ExpandMode="true";;
        H) echo "$usage"; exit;;
    esac
done

echo $ExpandMode

HG=hg19

### input information
INDELS=`readlink -f $INDELS` 
BAM=`readlink -f $BAM` 
BamNam=$(echo $INDELS | sed s/.txt// )
DIR=$BamNam'_snapshot' ## figure output folder
mkdir -p $DIR

###=========================================================
#ADDRESS="/home/local/ARCS/hq2130/src/IGV_2.3.68/"
ADDRESS="/home/local/ARCS/hq2130/src/IGV_2.3.68/"
IGVR=$ADDRESS"igv.sh" ## igv 
SCRF=$INDELS.txt
rm -rf $SCRF
touch $SCRF
printf "#! /bin/bash\n" >> $SCRF


kk=1
BAMS0="OO"
while read line
do
	## Step 1: write the tmp IGV igv_batch_script.txt
	
	NAME=`echo $line | cut -d ' ' -f 1`
	VALUE=`echo $line | cut -d ' ' -f 2`
	SAMPLE=`echo $line | cut -d  ' ' -f 3`
	SAMS=`echo $SAMPLE|tr "," "\n"`

	i=1
	for ONE in $SAMS;
 	do	
 		ONEBAM=`cat $BAM| grep $ONE$ ` # `grep $ONE "$BAM" | cut -f1`
 		#echo $ONE, $BAM, $ONEBAM
		if [[ "$ONEBAM" != "" ]]
		then
			if [[ $i -gt 1 ]]; then
				BAMS=`echo $BAMS,$ONEBAM`
			else
				BAMS=`echo $ONEBAM`
			fi
		fi
		let i+=1	
	done

	
	if [[  "$BAMS0" != "$BAMS" ]];then
		printf "new\n" >> $SCRF
		#printf "genome hg19\n"  >> $SCRF
		printf "genome $HG\n"  >> $SCRF
		printf "load  $BAMS\n" >> $SCRF
		printf "snapshotDirectory $DIR \n" >>  $SCRF
	fi
	
	BAMS0=$BAMS
	

	let "START=$VALUE - 50"
	let "END=$VALUE + 50"
	printf "goto chr$NAME:$START-$END \n" >> $SCRF
	printf "sort base \n" >> $SCRF
	if [[ $ExpandMode == "true" ]]; then
		#trackname=`echo $SAMPLE | cut -d ',' -f1`
		#printf "expand $trackname\n" >> $SCRF
		printf "expand \n" >> $SCRF

	else
		#printf "collapse \n" >> $SCRF
		printf "collapse \n" >> $SCRF
		

	fi


	printf "maxPanelHeight 400 \n" >> $SCRF
	#printf "snapshot $SAMPLE.$NAME.$VALUE.png \n" >> $SCRF
	SAMPLENAME=`basename $SAMPLE|grep -P -o '[a-zA-Z0-9-]+'`
	SAMPLENAME=`echo $SAMS|cut -d ' ' -f1 `
	#echo $SAMPLENAME
	printf "snapshot $SAMPLENAME.$NAME.$VALUE.png \n" >> $SCRF
	printf "\n" >> $SCRF
	
	let kk+=1
done < "$INDELS"

printf "exit \n" >> $SCRF
## run IGV
$IGVR  -g $HG -b $SCRF	
#rm $SCRF	
	
