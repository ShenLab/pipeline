#!/usr/bin/env python

# The purpose of this script is to filter a VCF file to remove all variants that have no alternate allele calls. Primarily it is for use with vcfs that have be subsetted by sample
#	-v/--vcf	<required>	The script requires an input vcf
#	-o/--out	<required>	The use should specify a base name for the output files. The script outputs the filtered results a vcf file. The script also produces a log file.

import os
from optparse import OptionParser
parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-o", "--out", dest="OutputFileName",help="user specified name of file for output to be written", metavar="OutputFileName")
(options, args) = parser.parse_args()

#Assign input and output files
VCF=open(options.VCFfile,'r')
BaseName=str(options.OutputFileName)
VcfOutputFilename=BaseName+'.vcf'
LogOutputFilename=BaseName+'.log'
Outvcf=open(VcfOutputFilename,'w')
Outlog=open(LogOutputFilename,'w')

#start log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("VCF filtering log: "+TimeNow+"\n")
Outlog.write("VCF: "+"\n")


# Read VCF file
BadGeno=['0/0', './.']
OrigCount=0
FiltCount=0
for line in VCF:
	OrigCount=OrigCount+1
	if '#' in line:
		Outvcf.write(line)
	if line.startswith("#CHROM"):
		testsamp=line.split('\t')
		SamNum=len(testsamp)
	if '#' not in line:
		# Variant must first pass 1KG and GO-ESP frequencies, MQ0 threshold, and be exonic
		linelist=line.split("\t")
		
		# Get Genotypes
		GenotypeString=[ linelist[i].strip() for i in range(9,SamNum) ]
		GenotypeList=[ i.split(':') for i in GenotypeString ]
		Genotypes=[ GenotypeList[i][0] for i in range(0,len(GenotypeString)) ]
		# Check if non-ref allele is present
		if not all(i in BadGeno for i in Genotypes):
			Outvcf.write(line)
			FiltCount=FiltCount+1
Outlog.write("Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("Number of variants left: "+str(FiltCount)+"\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Filtering Finished: "+TimeNow+"\n")
print "Filtering finished "+TimeNow
