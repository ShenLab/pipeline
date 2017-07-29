#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Annotate the variant list
# Variant list according to PSAP report.
#=========================================================================

import time
import argparse
import yaml
import csv
import os
import pprint
from MyVCF import *
import pandas as pd 

#HEADERS_PSAP_TO_BE_POP = ['Chr','Start','Ref','Alt','Gene.wgEncodeGencodeBasicV19', 'Func.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19','AAChange.wgEncodeGencodeBasicV19','mac63kFreq_ALL','1000g2014sep_all','esp6500si_all']
HEADERS_PSAP_TO_BE_POP = ['Gene.wgEncodeGencodeBasicV19', 'Func.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19','AAChange.wgEncodeGencodeBasicV19','mac63kFreq_ALL','1000g2014sep_all','esp6500si_all']
HeaderFromVCF_P1 = ['Chrom', 'Pos', 'Ref', 'Alt']
CSV_HEADER = ['Gene', 'GeneName', 'Allele Count', 'GeneFunc', 'ExonicFunc', 'AAchange', 'ExAC_ALL', 'gnomAD_genome_ALL', '1KG', 'VarType', 'MetaSVM', 'CADD13', 'PP2', 'MCAP', 'mis_z', 'lof_z', 'pLI', 'pRec','HeartRank', 'LungRank', 'BrainRank', 'Mappability', 'Filter']

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--report', type=str, help='PSAP Report')
	parser.add_argument('-v', '--vcf', type=str, help='VCF record')
	parser.add_argument('-p', '--ped', type=str, help='Pedigree')
	parser.add_argument('-o', '--out', type=str, help='Output Prifix')
	args = parser.parse_args()
	if args.out == None:
		args.out = args.ped.rstrip('ped')
	return args.report, args.vcf, args.ped, args.out

class PSAP_REPORT(object):
	"""docstring for PSAP_REPORT"""
	def __init__(self, PsapReportTxT, VarDict, Genes, OUT):
		self.VarDict = VarDict
		self.Genes = Genes.Genes
		self.fin = open(PsapReportTxT, 'rb')
		self.FullHeader = self.fin.readline().strip().split('\t')
		TmpHeader = self.FullHeader
		print TmpHeader
		pop_idxes = [TmpHeader.index(item) for item in HEADERS_PSAP_TO_BE_POP]
		pop_idxes.sort(reverse=True)
		print pop_idxes
		self.pop_idxes = pop_idxes
		self.Pop_non_display(TmpHeader)
		print TmpHeader
		self.Header = HeaderFromVCF_P1 + TmpHeader[4:]
		self.Header += CSV_HEADER
		self.Header += [x+'.GT' for x in self.VarDict.values()[0].headers[9:]]

		self.OutFil = open(OUT+'PSAP.csv','wb')
		self.Writer = csv.writer(self.OutFil)

	def RunReport(self):
		self.Writer.writerow(self.Header)
		for l in self.fin:
			row = l.strip().split('\t')
			self.Pop_non_display(row)
			record = Record()
			record.read_PsapReport(row, self.Header)
			record.fetch_VCF(self.VarDict)
			P1 = [record.VCF.Chrom, record.VCF.Pos, record.VCF.Ref, record.VCF.Alt]
			others = record.FormRecord(self.Genes)
			res = P1 + row[4:] + others
			self.Writer.writerow(res)
		self.OutFil.close()

	def Pop_non_display(self, row):
		for idx in self.pop_idxes:
			row.pop(idx)

class Record():
	def __init__(self):
		pass

	def read_PsapReport(self, row, header):
		self.CHROM = row[0]
		self.POS = row[1]
		self.REF = row[2]
		self.ALT = row[3]


	def fetch_VCF(self, vcf_dict):
		try:
			key = '{}:{}:{}:{}'.format(self.CHROM, self.POS, self.REF, self.ALT)
			#if self.ALT != '-':
			#	self.VCF = vcf_dict[self.CHROM + ':' + self.POS]
			#elif self.ALT == '-':
			#	self.VCF = vcf_dict[self.CHROM + ':' + str(int(self.POS) - len(self.REF) + 1)]
			self.VCF = vcf_dict[key]
		except KeyError:
			print "KeyError while fetch record in VCF, {}:{}:{}:{}".format(self.CHROM, self.POS, self.REF, self.ALT)
			exit()

	def FormRecord(self, genesocre):
		Gene = self.VCF.Info['Gene.refGene'][0]
		Gene = FixGene(Gene)
		try:
			GeneName = genesocre[Gene].Name
		except:
			GeneName = '.'
			tmp = GENE(Gene)
			genesocre[tmp.Symbol] = tmp
		try:
			AC = ','.join(self.VCF.Info['AC'])
		except:
			AC = ','.join(self.VCF.Info['MLEAC'])
		GeneFunc = ','.join(self.VCF.Info['Func.refGene'])
		ExonicFunc = ','.join(self.VCF.Info['ExonicFunc.refGene'])
		#AAchange = ','.join(self.VCF.Info['AAChange'])
		AAchange = ','.join(self.VCF.Info['AAChange.refGene'])
		ExAC_ALL = ','.join(str(AF(x)) for x in self.VCF.Info['ExAC_ALL'])
		gnomAD_genome_ALL = ','.join(str(AF(x)) for x in self.VCF.Info['gnomAD_genome_ALL'])
		MCAP = ','.join(self.VCF.Info['MCAP'])
		MetaSVM = ','.join(self.VCF.Info['MetaSVM_pred'])
		REVEL = ','.join(self.VCF.Info['revel'])
		CADD = ','.join(self.VCF.Info['CADD_phred'])
		PP2 = ','.join(self.VCF.Info['Polyphen2_HDIV_pred'])
		VarType = []
		for a,b,c,d,e in zip(GeneFunc.split(','), ExonicFunc.split(','), MetaSVM.split(','), CADD.split(','), PP2.split(',')):
			VarType.append(self.VCF.GetVarType(a,b,c,d,e))
		VarType = ','.join(VarType)
		_1KG = ','.join(str(AF(x)) for x in self.VCF.Info['1000g2015aug_all'])
		mis_z = genesocre[Gene].Mis_z
		lof_z = genesocre[Gene].Lof_z
		pLI = genesocre[Gene].pLI
		pRec = genesocre[Gene].pRec
		HeartRank = genesocre[Gene].DiaphragmRank
		LungRank = genesocre[Gene].LungRank
		BrainRank = genesocre[Gene].MouseBrainRank

		#return '\t'.join([self.CHROM, self.POS, self.REF, self.ALT, self.Gene, self.GeneName, self.GeneFunc, self.ExonicFunc, self.AAchange, self.ExACfreq, self.gnomAD_genome, self._1KGfreq, self.MetaSVM, self.CADD, self.PolyPhen2, self.pLI, self.misZ, self.lofZ, self.pREC, self.HRank, self.LRank, self.BRank, self.Proband, self.Father, self.Mother, self.DzModel, self.Pval, self.Flag])
		INFO = [Gene, GeneName, AC, GeneFunc, ExonicFunc, AAchange, ExAC_ALL, gnomAD_genome_ALL, _1KG, VarType, MetaSVM, CADD, PP2, MCAP, mis_z, lof_z, pLI, pRec, HeartRank, LungRank, BrainRank]
		GTs = self.VCF.getSampleGenotypes()
		return INFO + [self.VCF.Filter] + GTs

def FixGene(GeneName):
	if GeneName == 'BIVM-ERCC5-ERCC5':
		return 'ERCC5'
	else:
		return GeneName

def GenerateReport(psap_report, vcf, ped, OUT):
	f_psap = open(psap_report, 'rb')
	f_vcf = open(vcf, 'rb')
	#Ped = Pedigree(ped)
	Ped = None
	genescore = GENE_ANNOTATION()
	vcf_dict = GetVCF(f_vcf, Ped, genescore, OUT)
	Report = PSAP_REPORT(psap_report, vcf_dict, genescore, OUT)
	Report.RunReport()

def GetVCF(f_vcf, Ped, gene_score, OUT):
	print 'Reading VCF'
	stime = time.time()
	res = {}
	#INDEL_OUT = open(OUT + 'INDEL.csv', 'wb')
	#Writer = csv.writer(INDEL_OUT, delimiter=',')
	#Writer.writerow(CSV_HEADER)
	#Filters = Parse_YAML(ALL_FILTER)
	for l in f_vcf:
		if l.startswith('##'):
			continue
		elif l.startswith('#'):
			vcf_header = l.strip().split('\t')
		else:
			var = Variant(l, vcf_header)
			Add_var(var, res)
			#_pass, Proband, Father, Mother= var.ProbandisINDEL(Ped, Filters)
			#if _pass:
			#	Writer.writerow(var.OutAsCSV(Proband, Father, Mother, Ped, gene_score))

	print 'Finished Reading VCF %.3f' % -(stime - time.time())
	return res

def Add_var(var, res):
	Chrom = var.Chrom
	Pos = var.Pos
	Ref = var.Ref
	Alts = var.Alts
	for Alt in Alts:
		key = '{}:{}:{}:{}'.format(Chrom, Pos, Ref, Alt)
		res[key] = var
		pos, ref, alt = TrimAllele(Pos, Ref, Alt)
		if pos != None:	
			key = '{}:{}:{}:{}'.format(Chrom, pos, ref, alt)
			res[key] = var

def TrimAllele(Pos, Ref, Alt):
	# Left Trim
	pos, ref, alt = int(Pos), Ref, Alt
	while 1:
		if ref == '':
			ref = '-'
			pos -= 1
			break
		if alt == '' or alt == None:
			alt = "-"
			break
		if ref[0] == alt[0]:
			ref = ref[1:]
			alt = alt[1:]
			pos += 1
		else:
			return None,None,None	
	#if not (len(ref) == len(alt) and len(ref) == 1):
	#	print Pos, Ref, Alt
	#	print pos, ref, alt
	return str(pos), ref, alt

def main():
	PSAP, VCF, PED, OUT = GetOptions()
	GenerateReport(PSAP, VCF, PED, OUT)
	print 'Done'
	return


if __name__ == '__main__':
	main()
