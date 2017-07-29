#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Filters on VCF file
# Load an yml file and perform certein filters.
#=========================================================================

import time
import gzip
import argparse
import yaml
import csv
import re
import os
import pprint

PGLINE="##ExmFilter.py"
VQSR = re.compile('([\d.]+)')
pp = pprint.PrettyPrinter(indent=4)

class ExmFilter:
	def __init__(self, VCFName, YMLFilters, Ped=None, OutName=None, Debug=False):
		self.VCFName = VCFName
		if OutName == None:
			self.OutName = re.sub('\.gz$', '', self.VCFName)
			self.OutName = re.sub('\.vcf$', '', self.OutName)
			self.OutName = self.OutName + '.RareCoding.vcf'
		else:
			self.OutName = OutName
		print self.OutName
		self.Filter = YMLFilters
		self.Pedigree = Ped
		self.isTrio = self.Pedigree.isTrio()
		self.Debug = Debug
	def run(self):
		fin = GetVCF(self.VCFName)
		fout = open(self.OutName, 'wb')
		if self.Debug:
			ferr = open('FilteredOut.'+self.OutName, 'wb')
		for l in fin:
			if l.startswith('##'):
				fout.write(l)
				if self.Debug:
					ferr.write(l)
			elif l.startswith('#'):
				fout.write(PGLINE)
				fout.write(l)
				if self.Debug:
					ferr.write(PGLINE)
					ferr.write(l)
				self.Header = l.strip().split('\t')
			else:
				var = VARIANT(l)
				FLAG =  var.FilterByCriteria(self.Filter, self.Pedigree, self.Header)
				if FLAG == True:
					fout.write(l)
				elif self.Debug:
					ferr.write(var.DebugOut(FLAG))
		fin.close()
		fout.close()
		if self.Debug:
			ferr.close()

class YML_Filter():
	def __init__(self, YML_FILTER_FILE):
		global PGLINE
		print YML_FILTER_FILE
		PGLINE = PGLINE + " YMLFilters: " +YML_FILTER_FILE + "\n"
		with open(YML_FILTER_FILE, 'rb') as ymlfile:
			self.cfg = yaml.load(ymlfile)
		self.INFO = self.cfg.get('INFO', dict([]) )
		if self.INFO == None:
			self.INFO = dict([])
		self.READS = self.cfg.get('READS', dict([]) )
		if self.READS == None:
			self.READS = dict([])
		self.Filter = self.cfg.get('FILTER', dict([]) )
		if self.Filter == None:
			self.Filter = dict([])
		self.SNP = self.cfg.get('SNP', dict([]) )
		if self.SNP == None:
			self.SNP = dict([])
		self.INDEL = self.cfg.get('INDEL', dict([]) )
		if self.INDEL == None:
			self.INDEL = dict([])

	def show(self):
		pp.pprint(self.INFO)
		pp.pprint(self.READS)
		pp.pprint(self.FILTER)
		pp.pprint(self.SNP)
		pp.pprint(self.INDEL)

class INDIVIDUAL:
	def __init__(self, line):
		self.llist = line.strip().split('\t')
		self.FamID, self.SampleID, self.FatherID, self.MotherID, self.Gender, self.Affected = self.llist[0:6]

class PEDIGREE():
	def __init__(self, PedFil):
		print PedFil
		self.fin = open(PedFil,'rb')
		self.Header = self.fin.readline().strip().split('\t')
		self.individuals = []
		for l in self.fin:
			indi = INDIVIDUAL(l)
			self.individuals.append(indi)
	def isTrio(self):
		self.Proband, self.Father, self.Mother = None, None, None
		# Search for Proband.
		try:
			for indi in self.individuals:
				if indi.llist[self.Header.index('Relationship')] == 'Proband':
					self.Proband = indi.SampleID
					self.Father = indi.FatherID
					self.Mother = indi.MotherID
					if self.Father != '0' and self.Mother != '0':
						return True
					else:
						return False
		except ValueError:
			print "Proband not int header, Try to infer Proband by Affected"
			Max = 0
			Max_Candidate = None
			Tmp = 0
			if len(self.individuals) == 1: # Only one individual, this is proband
				self.Proband = self.individuals[0].SampleID
				print "Proband should be", self.Proband
				return
			for indi in self.individuals:
				if indi.Affected == '2': #Maybe a Proband
					if indi.FatherID != '0':
						Tmp += 1
					if indi.MotherID != '0':
						Tmp += 1
					if Tmp == 2: # If have father, mother and is affected, assume this is proband.
						self.Proband = indi.SampleID
						self.Father = indi.FatherID
						self.Mother = indi.MotherID
						print "Proband should be", self.Proband
						return True
					else: 
						if Tmp > Max:
							Max = Tmp
							Max_Candidate = indi
			if Max_Candidate != None:
				self.Proband = Max_Candidate.SampleID
				self.Father = Max_Candidate.FatherID
				self.Mother = Max_Candidate.MotherID
			else:
				self.Proband = None
				self.Father = None
				self.Mother = None
				print "This Pedigree don't have Proband."
				exit()
			print "Proband should be", self.Proband
			return False

class GENOTYPE():
	# GT:AD:DP:GQ:PL    0/0:7,0:7:18:0,18,270
	def __init__(self, Format, Genotype):
		self.Format = Format.split(':')
		self.Genotype = Genotype.split(':')
		if self.Genotype[0] == './.':
			self.isCalled = False
			return
		self.isCalled = True
		self.GT = map(int, re.findall('[\d.]', self.Genotype[0]))
		self.Dict = {}
		for i, item in enumerate(self.Format):
			self.Dict[item] = self.Genotype[i]

class VARIANT():
	def __init__(self, record):
		record = record.strip().split('\t')
		self.Chrom, self.Pos, self.Id, self.Ref, self.Alt, self.Qual, self.Filter, self.Info_str, self.Format = record[0:9]
		self.Genotypes = record[9:]
		self.Alts = self.Alt.split(',')
		self.Alleles = [self.Ref] + self.Alts
		self.GetInfo()

	def GetInfo(self):
		self.Info = {}
		tmp = self.Info_str.split(';')
		for kv in tmp:
			try:
				k, v = kv.split('=')
				if k not in self.Info:
					self.Info[k] = v.split(',')
				else:
					self.Info[k].extend(v.split(','))
			except:
				pass

	def CheckGT(self, Genotype, Filters):
		if not Genotype.isCalled:
			return "NotCall"
		GT1, GT2 = Genotype.GT
		DP = int(Genotype.Dict['DP'])
		ADs =  Genotype.Dict['AD'].split(',')
		AD1, AD2 = ADs[GT1], ADs[GT2]
		PL_NotRef = Genotype.Dict['PL'].split(',')[0]
		if Filters.READS.get('min_proband_AD',None) != None:
			if (int(AD1) < int(Filters.READS['min_proband_AD'])) or (int(AD2) < int(Filters.READS['min_proband_AD'])):
				return "min_proband_AD"
		if Filters.READS.get('min_proband_alt_freq',None) != None:
			Alt_freq = float(AD2)/(DP)
			if (Alt_freq < float(Filters.READS['min_proband_alt_freq'])):
				return "min_proband_alt_freq"
		if Filters.READS.get('min_proband_PL', None) != None:
			if int(PL_NotRef) < Filters.READS['min_proband_PL']:
				return "min_proband_PL"
		return True

	def CheckInfo(self, idx, Filters):
		# =============================================================================
		# Filter on Exon
		if Filters.INFO.get('exon', None)  == True:
			if Filters.INFO['exon_flag'] != None:
				VarFunc = self.Info.get('Func.refGene',[None]*(len(self.Alts)+1))[idx]
				if VarFunc not in Filters.INFO['exon_flag']:
					return 'exon_flag'
			if Filters.INFO.get('excluded_ExonicFunc', None) != None:
				VarClass = self.Info.get('ExonicFunc.refGene',[None]*(len(self.Alts)+1))[idx]
				if VarClass in Filters.INFO['excluded_ExonicFunc']:
					return 'excluded_ExonicFunc'
		# =============================================================================
		# =============================================================================
		# Filter on AC
		if Filters.INFO.get('max_AC', None) != None:
			if int(self.Info.get('AC',[0]*(len(self.Alts)+1))[idx]) > int(Filters.INFO['max_AC']):
				return 'max_AC'
			if int(self.Info.get('MLEAC',[0]*(len(self.Alts)+1))[idx]) > int(Filters.INFO['max_AC']):
				return 'max_AC'
		# =============================================================================
		# =============================================================================
		# Filter on Population MAF        
		if Filters.INFO.get('max_ExAC', None) != None:
			if self.GetAF(self.Info['ExAC_ALL'][idx]) > Filters.INFO['max_ExAC']:
				return 'max_ExAC'
		if Filters.INFO.get('max_gnomAD', None) != None:
			if self.GetAF(self.Info['gnomAD_genome_ALL'][idx]) > Filters.INFO['max_gnomAD']:
				return 'max_gnomAD'
		if Filters.INFO.get('max_1KG', None) != None:
			if self.GetAF(self.Info['1000g2015aug_all'][idx]) > Filters.INFO['max_1KG']:
				return 'max_1KG'
		# =============================================================================
		# =============================================================================
		# Filter on Exclude Gene/Chrom 
		if Filters.INFO.get('excluded_gene', None) != None:
			for gene in Filters.INFO['excluded_gene']:
				if gene in self.Info['Gene.refGene'][idx]:
					return 'excluded_gene'
		if Filters.INFO.get('excluded_chrom', None) != None:
			if self.Chrom in Filters.INFO['excluded_chrom']:
				return 'excluded_chrom'
		# =============================================================================
		# Filter on GenomeRegion/Mappability        
		if Filters.INFO.get('max_seqdup', None) != None:
			segdupScore = self.Info.get('genomicSuperDups',[0])[0]
			if segdupScore != '.':
				segdupScore = re.search(
						'Score:(\d+?\.?\d+)', segdupScore).group(1)
				if float(segdupScore) >= Filters.INFO['max_seqdup']:
					return 'max_seqdup'
		if Filters.INFO.get('Mappability', None) != None:
			if float(self.Info.get('Mappability', [1])[0]) < Filters.INFO['Mappability']:
				return 'Mappability'
		return True

	def isSNP(self, idxAllele):
		if len(self.Ref) == 1 and len(self.Alleles[idxAllele + 1]) == 1:
			return True
		else:
			return False

	def isINDEL(self, idxAllele):
		return not self.isSNP(idxAllele)

	def CheckSNP(self, Filters):
		if Filters.SNP.get('max_FS', None) != None:
			if float(self.Info.get('FS', [0])[0]) > Filters.SNP['max_FS']:
				return 'max_FS'
		if Filters.SNP.get('min_QD',None) != None:
			if float(self.Info.get('QD', [100])[0]) < Filters.SNP['min_QD']:
				return 'min_QD'
		if Filters.SNP.get('min_ReadPosRankSum', None) != None:
			if float(self.Info.get('ReadPosRankSum',[100])[0]) < Filters.SNP['min_ReadPosRankSum']:
				return 'min_ReadPosRankSum'
		return True

	def CheckINDEL(self, Filters):
		if Filters.SNP.get('max_FS', None) != None:
			if float(self.Info.get('FS', [0])[0]) > Filters.SNP['max_FS']:
				return 'max_FS'
		if Filters.SNP.get('min_QD',None) != None:
			if float(self.Info.get('QD', [100])[0]) < Filters.SNP['min_QD']:
				return 'min_QD'
		if Filters.SNP.get('min_ReadPosRankSum', None) != None:
			if float(self.Info.get('ReadPosRankSum',[100])[0]) < Filters.SNP['min_ReadPosRankSum']:
				return 'min_ReadPosRankSum'
		return True

	def CheckFilter(self, Filters, idxAllele):
		if self.Filter == 'PASS' or self.Filter == '.':
			return True
		try:
			v1, v2 = VQSR.findall(self.Filter)
			if Filters.Filter.get('VQSRSNP', None) != None and self.isSNP(idxAllele):
				if float(v2) > Filters.Filter.get('VQSRSNP', None):
					return 'VQSR'
			elif Filters.Filter.get('VQSRINDEL', None) != None:
				if float(v2) > Filters.Filter.get('VQSRINDEL', None):
					return 'VQSR'
			return True
		except:
			print self.Filter

	# Find the Allele Index
	# Need: Pedigree Proband, 
	def FindAllele(self, Pedigree, Header):
		if Pedigree.Proband != None:
			idx_Genotype = Header.index(Pedigree.Proband) - 9
			Genotype = GENOTYPE(self.Format, self.Genotypes[idx_Genotype])
			if Genotype.isCalled:
				return int(Genotype.GT[1]), Genotype
			else:
				return None, Genotype
		else:
			print "No Proband Find, use First one as ref" 
			Genotype = GENOTYPE(self.Format, self.Genotypes[1])
			return 1, Genotype

	def FilterByCriteria(self, Filters, Pedigree, Header):
		idxAllele, Genotype = self.FindAllele(Pedigree, Header)
		if idxAllele == None:
			return 'NotCalled'
		else:
			idxAllele -= 1
		FLAG_INFO = self.CheckInfo(idxAllele, Filters)
		if FLAG_INFO != True:
			return FLAG_INFO
		FLAG_GT = self.CheckGT(Genotype, Filters)
		if FLAG_GT != True:
			return FLAG_GT
		FLAG_FILTER = self.CheckFilter(Filters, idxAllele)
		if FLAG_FILTER != True:
			return FLAG_FILTER
		if self.isINDEL(idxAllele):
			FLAG_INDEL = self.CheckINDEL(Filters)
			if FLAG_INDEL != True:
				return FLAG_INDEL
		else:
			FLAG_SNP = self.CheckSNP(Filters)
			if FLAG_SNP != True:
				return FLAG_SNP            
		return True

	def GetAF(self, Freq):
		if Freq == '.':
			return 0
		else:
			return float(Freq)

	def DebugOut(self, FLAG):
		if self.Filter == '.':
			Filter = FLAG
		else:
			Filter = '{},{}'.format(self.Filter, FLAG)
		return '\t'.join([self.Chrom, self.Pos, self.Id, self.Ref, self.Alt, self.Qual, Filter, self.Info_str, self.Format, '\t'.join(self.Genotypes)]) + '\n'

def GetVCF(VCFin):
	if VCFin.endswith('.gz'):
		return gzip.open(VCFin, 'rb')
	else:
		return open(VCFin, 'rb')

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', '--vcf', required=True, type=str, help='VCF record')
	parser.add_argument('-o', '--out', type=str, help='Output Prifix')
	parser.add_argument('-p', '--ped', type=str, help='Pedigree file of the VCF samples')
	parser.add_argument('-f', '--filter', default='/home/local/users/jw/CUMC/Exome-Filters-Jiayao/ALL_FILTER.yml' ,type=str, help='YML file for filter criteria')
	parser.add_argument('--debug', default=False, type=bool, help='Whether Output Filtered Variants for Debug Purpose')
	args = parser.parse_args()
	return args.vcf, args.out, args.filter, args.ped, args.debug

def main():
	VCFin, VCFout, Filters, PedFil, Debug = GetOptions()
	YMLFilters = YML_Filter(Filters)
	Ped = PEDIGREE(PedFil)
	instance = ExmFilter(VCFin, YMLFilters, Ped=Ped, OutName=VCFout, Debug=Debug)
	instance.run()
	print "Done"

if __name__ == '__main__':
	main()
