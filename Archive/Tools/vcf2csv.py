#!/home/local/users/jw/bin/python2.7
# read a vcf file and write as tsv file
from utils import *
import re
import os
import shutil
import csv
import gzip

class Variant():
	def __init__(self, record, headers):
		self.record = record.strip()
		record = record.strip().split('\t')
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = record[0:9]
		self.Chrom = CHROM
		self.Pos = POS
		self.Id = ID
		self.Ref = REF
		self.Alt = ALT
		self.Qual = QUAL
		self.Filter = FILTER
		self.GetInfo(INFO)
		self.Format = FORMAT
		self.GetSample(record, headers)

	def GetInfo(self, INFO):
		info = {}
		tmp = INFO.split(';')
		for item in tmp:
			item = item.split('=')
			if len(item) == 2:
				if item[0] in info:
					info[item[0]] = info[item[0]] + ',' + item[1]
				else:
					info[item[0]] = item[1]
			elif len(item) == 1:
				info[item[0]] = None
			else:
				print "Error Parsing INFO:", item
		self.Info = info

	def GetSample(self, record, headers):
		self.Sample = {}
		for i, sample in enumerate(headers[9:]):
			self.Sample[sample] = record[i + 9]

	def show_samples(self):
		for k, v in self.Sample.items():
			print '%s:%v' % (k, v)


def get_header(line):
	return line.strip().split('\t')


def variant2tsv(variant, header, info_num, Ind):
	tmp = [variant.Chrom, variant.Pos, variant.Id,
			variant.Ref, variant.Alt, variant.Qual, variant.Filter]
	for field in header[7:info_num + 7]:
		try:
			tmp.append(variant.Info[field])
		except KeyError:
			# print field
			if field == 'ExACfreq':
				tmp.append('0')
			else:
				tmp.append('.')
	tmp.append(variant.Format)
	for field in header[info_num + 8:]:
		# print field
		if Ind != None:
			try:
				tmp.append(variant.Sample[field])
			except KeyError:
				print field
				tmp.append('.')
		else:
			try:
				if field in header:
					tmp.append(variant.Sample[field])
			except KeyError:
				print field
				tmp.append('.')
	# print '\t'.join(tmp)+'\n'
	for i, item in enumerate(tmp):
		if item == None:
			tmp[i] = '.'
	return tmp
#return '\t'.join(tmp) + '\n'


def format_info(INFOs):
	res = []
	info = re.compile('##INFO=<ID=([^,]+),')
	for item in INFOs:
		res.append(info.match(item).group(1))
	# print res
	return res


Custome_info = ['Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'MLEAC', 'MLEAF', 'gnomAD_genome_ALL', 'gnomAD_exome_ALL', '1000g2015aug_all', 'ExAC_ALL', 'MetaSVM_pred', 'CADD13_PHRED', 'MCAP', '', 'avsnp147', 'Annotation', 'GeneticModels', 'ModelScore', 'Compounds', ]
#Custome_info = ['VarClass','AAChange','ESPfreq','1KGfreq','ExACfreq','PP2.hvar.prd','MetaSVMprd','CADDphred','CADDInDelraw','CADDInDelphred','VarFunc','GeneName']


def format_header(header, Infos, Ind):
	res = header[:7]
	res.extend(Infos)
	if Ind == None:
		res.extend(header[8:])
	else:
		res.append(header[8])
		for item in header[9:]:
			if item in Ind:
				res.append(item)

	return res, len(Infos)


def GetOptions():
	parser = OptionParser()
	parser.add_option('-v', '--vcf', dest='VCF',
			metavar='VCF', help='Input VCF file name')
	parser.add_option('-o', '--outname', dest='OutName',
			metavar='OutName', help='Name of Output table')
	parser.add_option('-d', '--dir', dest='Dir',
			help='InputDir, if given, all vcf in the dir will be converted to tsv file')
	parser.add_option('-i', '--ind', dest='Ind',
			help='Sample Id in header to keep in resutls')
	(options, args) = parser.parse_args()
	if options.OutName == None and options.VCF != None:
		options.OutName = GetBaseName(options.VCF)
	return options.VCF, options.OutName, options.Dir, options.Ind


def VCF2TABLE(vcf, out, Ind):
	print vcf
	if Ind != None:
		fin = open(Ind, 'rb')
		Ind = [item.strip() for item in fin.readlines()]
	print Ind
	fin = GetVCF(vcf)
	INFOs = []
	header = []
	for line in fin:
		if line.startswith('##'):  # Mate Info
			if line.startswith('##INFO'):  # Info
				INFOs.append(line)
		elif line.startswith('#'):  # header line, next line is variant.
			# print line
			header = get_header(line)
			Infos = format_info(INFOs)
			#Infos = Custome_info
			new_header, info_num = format_header(header, Infos, Ind)
			res = []
			#fout.write('\t'.join(new_header) + '\n')
		else:
			variant = Variant(line, header)
			res.append(variant2tsv(variant, new_header, info_num, Ind))
			#fout.write(variant2tsv(variant, new_header, info_num, Ind))

	#Trim Individuals that don't have any variants.
	idx_last = len(new_header)
	idx_format = new_header.index('FORMAT')
	for j in xrange(idx_last, idx_format, -1):
		#print j
		FLAG = False
		for i in xrange(0, len(res), 1):
			#print i
			#print res[i]
			#print res[i][j-1]
			GT = res[i][j-1].split(':')[0]
			if GT != '0/0' and GT != './.':
				FLAG = True
				break
		if not FLAG:
			new_header.pop(j-1)
			for i in xrange(0, len(res)):
				res[i].pop(j-1)
	# Write the result
	fout = open(out + '.csv', 'wb')
	Writer = csv.writer(fout)
	Writer.writerow(new_header)
	for row in res:
		Writer.writerow(row)

def MultiVcf2Table(DIR, Ind):
	DIR = os.path.abspath(DIR)
	vcfs = get_files(DIR, '.vcf')

	for vcf in vcfs:
		outname = GetBaseName(vcf)
		VCF2TABLE(vcf, outname, Ind)


def main():
	VCF, OUT, DIR, Ind = GetOptions()
	if DIR == None:
		VCF2TABLE(VCF, OUT, Ind)
	else:
		MultiVcf2Table(DIR, Ind)

def GetVCF(VCFin):
    if VCFin.endswith('.gz'):
        return gzip.open(VCFin, 'rb')
    else:
        return open(VCFin, 'rb')

if __name__ == '__main__':
	main()
