#!/home/local/users/jw/bin/python2.7
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Filter and gruop variants by Genotype.
#=========================================================================

from optparse import OptionParser
from Variant import INFO, Sample, Variant


class Pedigree

# each pedigree has a fam,proband,father,mother,out_name


class Pedigree():
	def __init__(self, line):
		self.fam, self.proband, self.father, self.mother
		tmps = line.strip().split('\t')[0:4]
		self.out_name = '_'.join([tmps[1:4]) + '.vcf'
		self.fam= tmps[0]
		self.proband= tmps[1]
		self.father= tmps[2] if len(tmps[2]) > 1 else None
		self.mother= tmps[3] if len(tmps[3]) > 1 else None

# Return a list of pedigree obj
def GetPedigrees(ped_file):
	fin = open(ped_file, 'rb')
	res= []
	fin.readline()
	for line in fin:
		res.append(Pedigree(line))
	return res

# Autosome Recessive: 1/1 | 0/1 | 0/1
# Autosome Dominate:  0/1 or 1/1 | 0/1 or 1/1 | 0/1 or 1/1
# X-linked: At X chrom:
def isAR(GTp, GTf, GTm):
	if GTp[0] == GTp[1] and GTp[0] in GTf and GTp[0] in GTm:
		return True

def isAD_denovo()

def Filter_by_GT(vcf, ped):
	pedigrees= GetPedigrees(ped)
	# Open file handle for each family
	# ped_files = []
	# for pedigree in pedigrees:
	#	ped_files.append(pedigree.out_name,'wb')
	# Store each category into 2D list
	ped_vars= []
	for pedigree in pedigrees:
		ped_vars.append([[], [], []])  # AR-AD
	if vcf.startswith('.gz'):
		fin = gzip.open(vcf, 'rb')
		outname= vcf.strip('.vcf.gz') + '.txt'
	else:
		fin = open(vcf, 'rb')
		outname= vcf.strip('.vcf') + '.txt'
	fout = open(outname, 'wb')

	for l in fin:
		if l.startswith('##'):
			continue
		elif l.startswith('#'):
			header= l.strip().split('\t')
		else:  # Parse record
			var = Variant(line, header)
			for fam_i, fam in enumerate(pedigrees):
				try:
					if (var.Sample[fam.proband].GT)



def GetOptions():
	parser= OptionParser()
	parser.add_option('-v', '--vcf', dest='VCF', metavar='VCF',
	                  help='vcf file to process.')
	parser.add_option('-p', '--ped', dest='PED', metavar='PED',
	                  help='pedigree file needed.')
	parser.add_option('-d', '--dir', dest='DIR', metavar='DIR',
	                  help='directory to save output files.')
	(options, args) = parser.parse_args()

	return options.VCF, options.PED

def main():
	vcf, ped = GetOptions()
	Filter_by_GT(vcf, ped)


if __name__ == '__main__':
	main()
