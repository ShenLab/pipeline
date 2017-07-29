#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# FamAnalysis.py
# Given A Cohort Pedigree file and Cohort VCF file, Generate Scripts for Family analysis, such as PSAP 
# And De novo Calling.
#========================================================================================================

import argparse
import os
import subprocess
import gzip
PrepareVCF = '$HOME/CUMC/Tools/PSAP_Family.sh'
COLLECT = '$HOME/CUMC/Tools/PSAP_CollectResults.py'
PREPARE_IGV = '$HOME/CUMC/Tools/PSAP_PrepareIGVInput.py'

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v','--vcf', type=str, required=True, help = 'VCF File')
	parser.add_argument('-p','--ped', type=str, required=True, help = 'Pedigree File')
	parser.add_argument('-d','--dir', type=str, help = 'work dir to generate results')
	parser.add_argument('--parallel', type=int, help='Number of Parallel to go. ')
	args = parser.parse_args()
	if args.dir == None:
		args.dir="./" + args.ped.rstrip(".ped") + ".PSAP"
	#return args.vcf, args.ped, args.dir
	return args

class PedigreeIndv:
	def __init__(self,FamID, SampleID, FatherID, MotherID, Sex, Affected, Relationship):
		self.FamID = FamID
		self.SampleID = SampleID
		self.FatherID = FatherID
		self.MotherID = MotherID
		self.Sex = Sex
		self.Affected = Affected
		self.Relationship = Relationship

class Pedigree:
	def __init__(self):
		self.Individuals = []
		self.isTrio = False

class FamAnalysis:
	def __init__(self, args):
		self.VCFname = os.path.abspath(args.vcf)
		self.PEDname = args.ped
		self.RESdir = os.path.abspath(args.dir)
		tmp = self.PEDname.split('/')[-1].rstrip('ped') + self.VCFname.split('/')[-1].rstrip('.gz').rstrip('vcf')
		self.ArrFil = 'Arr.' + tmp + 'list'  
		self.CmdFil = 'Cmd.' + tmp + 'sh'
		self.foutArr = open(self.ArrFil,'wb')
		self.foutCmd = open(self.CmdFil,'wb')
		self.Nparallel = args.parallel
	def run(self):
		Samples = self.GetIndividualList()
		WorkDir = self.RESdir
		fin = open(self.PEDname, 'rb')
		NewFam = 'None'
		for l in fin:
			if l.startswith('#'):
				header = l
				continue
			llist = l.strip().split()
			FamID = llist[0]
			if FamID != NewFam:
				if NewFam != 'None':
					self.ProcessOneFam(WorkDir, Samples, NewFam, Fam, self.VCFname)
					NewFam = FamID
					Fam = [header]
					Fam.append(l)
				else:
					NewFam = FamID
					Fam = [header]
					Fam.append(l)
			else:
				Fam.append(l)
		self.ProcessOneFam(WorkDir, Samples, NewFam, Fam, self.VCFname)
		#CheckAndSeperate('{}.list'.format(NewFam), Mapper)		
		fin.close()
		self.WriteScript()

	def GetIndividualList(self):
		with gzip.open(self.VCFname) as fin:
			for l in fin:
				if l.startswith('##'):
					continue
				elif l.startswith('#'):
					print l
					print l.strip().split('\t')[9:]
					return l.strip().split('\t')[9:]

	def ProcessOneFam(self, WorkDir, SampleList, NewFam, Fam, VCF):
		#print Fam
		for item in self.printSample(Fam):
			if item not in SampleList:
				os.chdir(WorkDir)
				return
		print NewFam
		subdir = WorkDir + '/' + NewFam
		if not os.path.isdir(subdir):
			os.mkdir(subdir)
		os.chdir(subdir)
		fout = open('{}.ped'.format(NewFam), 'wb')	
		for item in Fam:
			fout.write(item)
		fout.close()
		fout = open('{}.list'.format(NewFam), 'wb')	
		for item in self.printSample(Fam):
			fout.write(item+'\n')
		fout.close()
		SampleList = '{}.list'.format(NewFam)
		if VCF != None:
			#cmd = 'nohup {} {} {} {} &'.format(PrepareVCF, VCF, SampleList, subdir)
			cmd = '{}/{}\n'.format(subdir,SampleList)
			#print cmd
			#process = subprocess.Popen(cmd, shell=True)
			self.foutArr.write(cmd)
		os.chdir(WorkDir)

	def printSample(self, Fam):
		for l in Fam:
			if l.startswith('#'):
				continue
			else:
				yield l.strip().split('\t')[1]

	def WriteScript(self):
		self.foutCmd.write('#/bin/bash\n\n')
		self.foutCmd.write('CMD={}\n'.format(PrepareVCF))
		self.foutCmd.write('VCF={}\n'.format(self.VCFname))
		self.foutCmd.write('InpFil={}\n'.format(os.getcwd() + '/' + self.ArrFil))
		self.foutCmd.write('COLLECT={}\n'.format(COLLECT))
		self.foutCmd.write("Num_Job=`wc -l $InpFil | awk '{print $1}'`\n\n")
		self.foutCmd.write('seq $Num_Job | parallel -j {} --ETA bash $CMD -i $InpFil -v $VCF -a {}'.format(self.Nparallel, "{}"))
		self.foutCmd.write('python $COLLECT -d {} \n'.format(self.RESdir))
		self.foutCmd.write('python $PREPAREIGV -d {} \n'.format(PREPARE_IGV))
		self.foutCmd.close()

def main():
	#VCF, PED, DIR = GetOptions()
	args = GetOptions()
	instance = FamAnalysis(args)
	instance.run()
	return

if __name__=='__main__':
	main()
