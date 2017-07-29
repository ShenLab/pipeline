#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# ModifySampleName.py
#========================================================================================================

import argparse
import re

N0 = re.compile('COL-CHUNG_([A-Za-z0-9-]+)_')
N1 = re.compile('_([A-Za-z]+)')
N2 = re.compile('([0-9]+)_')

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v','--vcf',required=True, type=str, help = 'VCF file to be modify sample names')
	parser.add_argument('-o','--out',required=True, type=str, help = 'Output name of modify sample names VCF')
	args = parser.parse_args()
	
	return args.vcf, args.out

def ModifySampleID(samples):
	for i in range(len(samples)):
		P = N0.search(samples[i]).group(1)
		samples[i] = P

def main1():
	VCF, OUT = GetOptions()
	fin = open(VCF, 'rb')
	fout = open(OUT+'.vcf', 'wb')
	for l in fin:
		if l.startswith('##'):
			fout.write(l)
		elif l.startswith('#'):
			fout.write('##ModifySampleNameWithPythonScript ModifySampleName.py\n')
			headers = l.strip().split('\t')
			Part1 = headers[:9]
			samples = headers[9:]
			#print Part1
			ModifySampleID(samples)
			#print samples
			tmp = Part1 + samples
			fout.write('\t'.join(tmp)+'\n')
		else:
			fout.write(l)
	return

def main2():
	fin1 = open('Cardiomyopathy.vcf','rb')
	fin2 = open('CardiomyopathyALL.ped','rb')
	fout = open('Individuals.Filtered.list','wb')
	for l in fin1:
		if l.startswith('##'):
			pass
		elif l.startswith('#'):
			headers = l.strip().split('\t')
			Part1 = headers[:9]
			VCFsamples = headers[9:]
	fin1.close()
	PEDsamples = []
	for l in fin2:
		if l.startswith('#'):
			pass
		else:
			llist = l.strip().split('\t')
			PEDsamples.append(llist[1])
			fout.write(llist[1]+'\n')
	fin2.close()
			
	# samples in VCF not PED:
	print "samples in VCF not PED"
	for sample in VCFsamples:
		if sample not in PEDsamples:
			print sample,
	print
	print "samples in PED not in VCF:"
	for sample in PEDsamples:
		if sample not in VCFsamples:
			print sample,
	print

def main():
	main1()
	#main2()

if __name__=='__main__':
	main()
