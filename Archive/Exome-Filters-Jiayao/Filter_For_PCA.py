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
	parser.add_option("-c", "--control", dest="Control", metavar="Control",
			help="Control data"
			)
	(options, args) = parser.parse_args()
	VCFin = options.VCF
	VCFout = options.OutVCF
	if options.Control == None:
		Control = "/home/local/users/jw/resources/AncestryPCA/resources/1KG_AJ_Domi_PCAcontrol.vcf.gz"
	else:
		Control = options.Control

	return VCFin, VCFout, Control

class Controls:
	def __init__(self,Control):
		self.OneKGFil = Control
		self.Out = open("1KG.vcf",'wb')
	def LoadVar(self):
		self.Variants = {}
		fin = gzip.open(self.OneKGFil, 'rb')
		for l in fin:
			if l.startswith("#"):
				self.Out.write(l)
				continue
			llist = l.strip().split("\t")
			#self.Variants["{}:{}".format(llist[0],llist[1])] = llist[0:9]
			self.Variants["{}:{}".format(llist[0],llist[1])] = l

def Filter(VCFin, VCFout, Control):
	_1KG = Controls(Control)
	_1KG.LoadVar()
	#_1KG = _1KG.Variants
	if VCFin.endswith('.gz'):
		hand = gzip.open(VCFin, 'rb')
	else:
		hand = open(VCFin, 'rb')
	if VCFout.endswith('.gz'):
		fout = gzip.open(VCFout, 'wb')
	else:
		fout = open(VCFout, 'wb')
	Count_All = 0
	Count_Pass = 0
	for l in hand:
		if l.startswith('##'):
			fout.write(l)
			continue
		elif l.startswith('#'):
			HeaderLine = l
			fout.write(l)
			continue
		llist = l.strip().split('\t')
		Count_All += 1
		key = "%s:%s"%(llist[0],llist[1])
		# Drop the variants that not present in 1KG 
		if key not in _1KG.Variants: # This Variant not present in Control Data
			continue
		else:
			#_1KG.pop(key)
			_1KG.Out.write(_1KG.Variants[key])
			#if F_ExAC_Coding(llist[7], 0.05):
			fout.write(l)
			Count_Pass += 1
		if Count_All % 10000 == 0:
			print "Read %d Variants" % Count_All
	print "Finish Reading All %d Variants. %d Variants Pass the Filter" % (Count_All, Count_Pass)
	exit()
	print "Starts add missing 1KG Genotypes."
	LenOfGT = len(HeaderLine.strip().split("\t")) - 9 
	for k,v in _1KG.items():
		fout.write("\t".join(v)+"\t"+"\t".join(["0/0"]*LenOfGT)+"\n")
	fout.close()

def F_ExAC_Coding(INFO, cutoff):
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

	for v in infodict['1000g2015aug_all'].split(','):
		try:
			if float(v) >= cutoff:
				return True
		except ValueError:
			return False 
	return False


def main():
	VCFin, VCFout, Control = GetOptions()
	Filter(VCFin, VCFout, Control)


if __name__ == '__main__':
	main()
