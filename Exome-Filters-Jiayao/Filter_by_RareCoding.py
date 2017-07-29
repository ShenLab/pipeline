#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Filters on VCF file
#=========================================================================

from optparse import OptionParser
import gzip
import os
import re

def GetOptions():
	parser = OptionParser()
	parser.add_option("-v", "--vcf", dest="VCF",
			help="Input VCF file", metavar="VCFfile")
	parser.add_option("-o", "--outvcf", dest="OutVCF",
			help="Name of Output VCF file", metavar="OutVCF", default="Filterd.vcf")
	parser.add_option("-d", "--debug", dest="debug", metavar="Filters",
			help="Filters apply on the variants. Splited by ','"
			)
	(options, args) = parser.parse_args()
	VCFin = options.VCF
	VCFout = options.OutVCF
	if options.debug != None:
		Debug = True
	else:
		Debug = False

	return VCFin, VCFout, Debug 

def Filter(VCFin, VCFout, Debug):
	if VCFin.endswith('.gz'):
		hand = gzip.open(VCFin, 'rb')
	else:
		hand = open(VCFin, 'rb')
	if VCFout.endswith('.gz'):
		fout = gzip.open(VCFout, 'wb')
	else:
		fout = open(VCFout, 'wb')
	if Debug:
		f_debug = open("Filtered.vcf", 'wb')
	Count_All = 0
	Count_Pass = 0
	for l in hand:
		if l.startswith('##'):
			fout.write(l)
			continue
		elif l.startswith('#'):
			fout.write('##FiltersByRareCoing="ExAC_ALL:1e-2. RefGene.Func:exonic,splicing,exonic-splicing"\n')
			fout.write(l)
			continue
		llist = l.strip().split('\t')
		Count_All += 1
		CodingRegion = ['exonic', 'splicing', 'exonic-splicing']
		Flag = F_ExAC_Coding(llist[7], 1e-2, CodingRegion)
		if Flag == True:
			fout.write(l)
			Count_Pass += 1
		elif Debug :
			f_debug.write(Flag)
			f_debug.write(l)

		if Count_All % 10000 == 0:
			print "Read %d Variants" % Count_All
	fout.close()
	print "Finish Reading All %d Variants. %d Variants Pass the Filter" % (Count_All, Count_Pass)

def F_ExAC_Coding(INFO, ExACfreq, Region):
	infolist = INFO.split(';')
	infodict = {}
	for kv in infolist:
		kv = kv.split('=')
		if len(kv) == 2:
			k, v = kv
			if k in infodict:
				infodict[k] = infodict[k] + ',' + v
			else:
				infodict[k] = v

	Pass_region = False
	#print Region, infodict['Func.refGene'].split(',')
#	for item, v in zip(Region, infodict['Func.refGene'].split(',')):
#		print "||{}||{}||".format(item, v)
#		if item == v:
#			print "Yeah"
#			Pass_region = True
#			break
	VarRegions = set(infodict['Func.refGene'].split(','))
	Region = set(Region)
	Pass_region = ( len(VarRegions.intersection(Region)) >= 1 )
	
	if Pass_region == False:
		return "Not Pass Coding Region"

	for v in infodict['ExAC_ALL'].split(','):
		try:
			if float(v) <= ExACfreq:
				return True
		except ValueError:
			return True
	return "Not Pass ExAC Freq"


def main():
	VCFin, VCFout, Debug = GetOptions()
	Filter(VCFin, VCFout, Debug)


if __name__ == '__main__':
	main()
