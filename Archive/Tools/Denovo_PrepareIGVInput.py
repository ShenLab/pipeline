#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# DNM_PrepareIGVInput.py
# Prepare bamout script and variant table for IGV shot
# Input: Dir contain subdirs of PSAP results, (PSAP.csv, .ped file is need)
#        TXT file contains (sampleID,BAM) tuple
# Output: List of (sampleID Interval1 Interval2 ... Inverbaln)
#         Scripts to run BAMOUT.
#         Scripts to run IGV
#=========================================================================

import argparse
import csv
import re
import os


SampleName = re.compile('([A-Za-z0-9-]+).bam')
ProjectName = re.compile('\w+')
EXM_REF = '$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh'
BAMOUT_CMD = '$HOME/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.15b.BamOut.sh'
IGV_GENERATE_CMD = '$HOME/CUMC/Exome-IGV-Jiayao/Generate_IGV_plots.sh'
RUN = 'nohup'

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', '--var', type=str, required=True,
			help='denovo csv file')
	parser.add_argument('-p', '--ped', type=str, required=True,
			help='ped file')
	parser.add_argument('-b', '--bam', type=str, required=True,
			help='File contains bam locations. ')
	parser.add_argument('-l1', '--list1', type=str,
			help='list of sample and corresponding intervals, For Bamout')
	parser.add_argument('-l2', '--list2', type=str,
			help='list of variants and corresponding bam, For IGV')
	parser.add_argument('-s1', '--script1', type=str,
			help='Script to run GATK Bamout. ')
	parser.add_argument('-s2', '--script2', type=str,
			help='Script to run IGV. ')
	parser.add_argument('--parallel', type=int,
			default=20, help='Number of Parallel to go. ')
	args = parser.parse_args()
	if args.list1 == None:
		args.list1 = 'Denovo.bamout.input'
	if args.list2 == None:
		args.list2 = 'Denovo.IGV.varlist'
	if args.script1 == None:
		args.script1 = 'run_bamout.sh'
	if args.script2 == None:
		args.script2 = 'run_IGV.sh'
	return args

class DNM_PrepareIGVInput:
	def __init__(self, args):
		self.VARFil = os.path.abspath(args.var)
		self.PEDFil = os.path.abspath(args.ped)
		self.Bams = BamLocation(args.bam)
		self.VarList1 = os.path.abspath(args.list1)
		self.VarList1_hand = open(self.VarList1, 'wb')
		self.VarList2 = os.path.abspath(args.list2)
		self.VarList2_hand = open(self.VarList2, 'wb')
		self.Script1 = open(args.script1 , 'wb')
		self.Script2 = open(args.script2 , 'wb')
		self.parallel = args.parallel
	
	def run(self):
		pedigree = Pedigree(self.PEDFil)
		self.variants = []
		with open(self.VARFil, 'rb') as csvfile:
			reader = csv.reader(csvfile)
			header = reader.next()
			for row in reader:
				var = Variant(header, row)
				var.addPed(pedigree)
				var.addBam(self.Bams)
				self.variants.append(var)
				self.VarList2_hand.write(var.OutVar())
		self.FlushBamout()
		self.FlushIGV()

	def ReduceVarian2Sample(self):
		res = {}
		for variant in self.variants:
			if variant.Sample not in res:
				res[variant.Sample] = Sample(variant)
			else:
				res[variant.Sample].AddVar(variant)
		print res
		print res.values()
		return res.values()

	def FlushBamout(self):
		# Write Bamout List	
		Samples = self.ReduceVarian2Sample()
		for sample in Samples:
			for line in sample.FlushBamout():
				self.VarList1_hand.write(line)
		# Write Bamout Script	
		fout = self.Script1
		fout.write('REF={}\n'.format(EXM_REF))
		fout.write('CMD={}\n'.format(BAMOUT_CMD))
		fout.write('InpFil={}\n'.format(self.VarList1))
		fout.write("NumJob=`wc -l $InpFil|cut -f1 -d ' '`\n\n")
		fout.write("mkdir -p Bamouts;cd Bamouts\n")
		fout.write(
				"seq $NumJob| parallel -j {} --eta $CMD -i $InpFil -a {} -r $REF\n\n".format(self.parallel, "{}"))
		fout.write("cd ../;find `pwd`/Bamouts/ -name \"*.bamout.bam\" > Bamout.bam.list\n")
		fout.write("cat {} {} > ALL.bam.list\n".format(self.Bams.bamlocationfile, "Bamout.bam.list"))
		fout.close()

	def FlushIGV(self):
		fout = self.Script2
		fout.write(
				'{} -b ALL.bam.list -v {}'.format(IGV_GENERATE_CMD, self.VarList2))
		fout.close()


class Sample:
	def __init__(self, variant):
		self.Sample = variant.Sample
		self.Father = variant.Father
		self.Mother = variant.Mother
		self.SampleBam = variant.SampleBam
		self.FatherBam = variant.FatherBam
		self.MotherBam = variant.MotherBam
		self.SampleBamout = variant.SampleBamout 
		self.FatherBamout = variant.FatherBamout
		self.MotherBamout = variant.MotherBamout
		# self.variants = []
		chrom, start, end = variant.Chrom, variant.start, variant.end 
		self.Intervals = ['{}:{}-{}'.format(str(chrom), str(start), str(end))]
	def AddVar(self, variant):
		# self.vairants.append(variant)
		chrom, start, end = variant.Chrom, variant.start, variant.end
		self.Intervals.append('{}:{}-{}'.format(str(chrom), str(start), str(end))) 
	def FlushBamout(self):
		# Format is
		# BAMPATH	-L	VAR1	-L VAR2	 .... -L VARN
		# SampleCMD =  ''.format(self.SampleBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
		SampleCMD =  '{}\t{} \n'.format(self.SampleBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
		FatherCMD =  '{}\t{} \n'.format(self.FatherBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
		MotherCMD =  '{}\t{} \n'.format(self.MotherBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
		return SampleCMD, FatherCMD, MotherCMD

class BAM:
	def __init__(self, fpath):
		self.FullPath = fpath.strip()
		self.BamName = self.FullPath.split('/')[-1]
		self.BamOutName = '.'.join(self.BamName.split('.')[:-1]) + ".bamout.bam"
		self.Sample = SampleName.search(self.BamName).group(1)

class BamLocation:
	def __init__(self, bamlocationfile):
		self.bamlocationfile = os.path.abspath(bamlocationfile)
		fin = open(bamlocationfile, 'rb')
		self.Bams = {}
		for l in fin:
			if 'bamout' in l:
				continue
			bam = BAM(l)
			self.Bams[bam.BamName] = bam
	def search(self,sample):
		for k,v in self.Bams.items():
			if sample in v.FullPath:
				return v
		return None

class Variant:
	def __init__(self, header, row):
		self.Chrom = row[header.index('Chrom')]
		self.Pos = int(row[header.index('Pos')])
		self.Sample = row[header.index('Sample')]
		self.start = max(0, self.Pos - 150)
		self.end = self.Pos + 150

	def addPed(self, pedigree):
		self.Father = pedigree.Probands[self.Sample].Father
		self.Mother = pedigree.Probands[self.Sample].Mother

	def addBam(self, bamlocation):
		#print self.Sample
		for bam in bamlocation.Bams.values():
			if self.Sample == bam.Sample:
				self.SampleBam = bam
				self.SampleBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
			if self.Father == bam.Sample:
				self.FatherBam = bam
				self.FatherBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
			if self.Mother == bam.Sample:
				self.MotherBam = bam
				self.MotherBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
		try:
			for var in [self.SampleBam, self.SampleBamout, self.FatherBam, self.FatherBamout, self.MotherBam, self.MotherBamout ]:
				if var not in locals():
					print 'Cant find {} in locals.'.format(var)
					exit()
		except:
			print self.Sample
	def OutVar(self):
		return '{}\t{}\t{}\n'.format(self.Chrom, self.Pos, ','.join([self.SampleBam.BamName, self.SampleBamout, self.FatherBam.BamName, self.FatherBamout, self.MotherBam.BamName, self.MotherBamout]))
	def Out2L(self):  # To run_Bamout.sh
		return "-L {}:{}-{}".format(self.Chrom, self.start, self.end)
	def Out2I(self):  # To run_IGV.sh
		return "{}\t{}".format(self.Chrom, self.Pos)

class Proband:
	def __init__(self, Sample, Father, Mother):
		self.Sample = Sample
		self.Father = Father
		self.Mother = Mother


class Pedigree:
	def __init__(self, pedfile):
		fin = open(pedfile, 'rb')
		self.Probands = {}
		for l in fin:
			if l.startswith('#'):
				continue
			llist = l.strip().split('\t')
			fam, sample, father, mother, sex, phenotype = llist[:6]
			if father != "0" and mother != "0" and phenotype == "2":  # Proband
				# print sample, father, mother
				self.Probands[sample] = Proband(sample, father, mother)

def main():
	args = GetOptions()
	instance = DNM_PrepareIGVInput(args)
	instance.run()
	return


if __name__ == '__main__':
	main()
