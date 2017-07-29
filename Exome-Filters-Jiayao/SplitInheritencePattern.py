#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# Split vcf file by inheritence pattern given by genmod
#========================================================================================================

from optparse import OptionParser
import utils
import re
import os
import shutil

def GetOptions():
	parser = OptionParser()
	parser.add_option('-d','--dir',dest = 'dir', metavar = 'dir', help = 'InputDir contains VCF file and pedigree files')
	parser.add_option('-v','--vcf',dest = 'vcf', help = 'Inhretence annotated VCF file')
	parser.add_option('-p','--ped',dest = 'ped', help = 'Pedigree file')
	
	(options,args) = parser.parse_args()
	
	return options.dir, options.vcf, options.ped

class Individual():
	def __init__(self, List):
		self.Fam, self.Sample, self.Father, self.Mother, self.Gender, self.Pheno = List[:6]

class Pedigree():
	def __init__(self, PedFil):
		fin = open(PedFil,'r')
		individuals = []
		for l in fin:
			if l.startswith('#'):
				continue
			indi = l.strip().split('\t')
			individuals.append(Individual(indi))
		self.Proband, self.Father, self.Mother = None, None, None
		for ind in individuals:
			if ind.Fam in ind.Sample:
				self.Proband = ind
		for ind in individuals:
			if self.Proband.Father == ind.Sample :
				self.Father = ind
			if self.Proband.Mother == ind.Sample :
				self.Mother = ind
class Format():
    #GT:AD:DP:GQ:PL 0/0:7,0:7:18:0,18,270
    def __init__(self,Format,tmp):
		self.formats = Format.split(':')
		self.info = tmp.split(':')
		self.GT = re.findall('[\d.]',self.info[0])
		self.fmt = {}
		for i,item in enumerate(self.formats):
			self.fmt[item] = self.info[i]
		if 'AD' in self.fmt:
			self.fmt['AD'] = self.fmt['AD'].split(',')

def MatchVcfPed(vcfs,peds):
	res = []
	for vcf in vcfs:
		for ped in peds:
			if utils.GetBaseName(vcf) == utils.GetBaseName(ped): 
				print utils.GetBaseName(vcf), utils.GetBaseName(ped)
				res.append((vcf,ped))
				break
	print "Finded vcf-ped pairs:"
	for v,p in res:
		print v,p
	print '\n'
	return res

def getGT(GTstring):
	return re.findall('[\d.]',GTstring.split(':')[0])

def CompoundTrios(CH_buffer, fout_CH, headerList, Ped):
	Flag1 = False
	Flag2 = False
	res = []
	Proband_idx = headerList.index(Ped.Proband.Sample)
	Father_idx = headerList.index(Ped.Father.Sample)
	Mother_idx = headerList.index(Ped.Mother.Sample)
	for l in CH_buffer:
		llist = l.strip().split('\t')
		GT_Proband = getGT(llist[Proband_idx])
		if GT_Proband[0] == GT_Proband[1]:
			continue
		GT_Father = getGT(llist[Father_idx])
		GT_Mother = getGT(llist[Mother_idx])
		if GT_Proband[1] in GT_Father and GT_Proband[1] not in GT_Mother:
			Flag1 = True
			res.append(l)
		if GT_Proband[1] in GT_Mother and GT_Proband[1] not in GT_Father:
			Flag2 = True
			res.append(l)
	if Flag1 and Flag2:
		for l in res:
			fout_CH.write(l)
def CompoundDual(CH_buffer,fout_CH, Ped):
	Flag1 = False
	Flag2 = False
	res = []
	for l in CH_buffer:
		llist = l.strip().split('\t')
		GT_Proband = getGT(llist[Proband_idx])
		if GT_Proband[0] == GT_Proband[1]:
			continue
		GT_Parent = getGT(llist[Parent_idx])
		if GT_Proband[1] in GT_Parent:
			Flag1 = True
			res.append(l)
		if GT_Proband[1] not in GT_Parent:
			Flag2 = True
			res.append(l)
	if Flag1 and Flag2:
		for l in res:
			fout_CH.write(l)
def CompoundQuart(CH_buffer,fout_CH, Ped):
	return

def LookCompoundHeter(CH_buffer,fout_CH,RecordLen,Format_idx, headerList, Ped):
	if RecordLen < 2:
		return
	if RecordLen - Format_idx - 1 == 1: # Singleton
		return
	if RecordLen - Format_idx - 1 == 2: # Dual
		Proband_idx, Parent_idx = [Format_idx+1, Format_idx+2]
		CompoundDual(CH_buffer,fout_CH,Proband_idx,Parent_idx)
	if RecordLen - Format_idx - 1 == 3: # Trios
		Proband_idx,Father_idx,Mother_idx = [Format_idx+i+1 for i in range(3)]
		CompoundTrios(CH_buffer,fout_CH, headerList, Ped)
	if RecordLen - Format_idx - 1 == 4: # quart	
		CompoundQuart()
		
	return

def selectAD_DN(llist, headerList, Ped, Model, Format_idx):
	Flag_pass = True
	if 'dn' in Model:
		try:
			Proband_idx = headerList.index(Ped.Proband.Sample)
			Father_idx = headerList.index(Ped.Father.Sample)
			Mother_idx = headerList.index(Ped.Mother.Sample)
			_format = llist[Format_idx]
			Proband = Format(_format, llist[Proband_idx])
			Father = Format(_format, llist[Father_idx])
			Mother = Format(_format, llist[Mother_idx])

			#print Proband.fmt['AD'], int(Proband.fmt['GT'][1]), float(Proband.fmt['AD'])
			if Proband.GT[1] == '.':
				return False
			ADDP = float(Proband.fmt['AD'][int(Proband.GT[1])]) / float(Proband.fmt['DP'])
			if (ADDP < 0.2) or (float(Proband.fmt['GQ']) < 30) or (Proband.fmt['AD'][int(Proband.GT[1])] < 3) or (int(Father.fmt['DP'])<10) or (int(Mother.fmt['DP'])<10):
				Flag_pass = False
				#print Proband.fmt['AD'],Proband.fmt['DP'], Proband.fmt['GQ']
				#print "Allele Balance: %.3f"%ADDP
			else:
				#print Proband.fmt['AD'],Proband.fmt['DP']
				#print "Allele Balance: %.3f"%ADDP
				pass	
		except Exception as e:
		#except IOError as e:
			print e, Proband.fmt
			Flag_pass = False
	return Flag_pass
	
def selectAR(llist, headerList, Ped, Model, Format_idx):
	Flag_pass = True
	try:
		Proband_idx = headerList.index(Ped.Proband.Sample)
		_format = llist[Format_idx]
		Proband = Format(_format, llist[Proband_idx])
		if Proband.GT[1] == '.':
			return False
	except Exception as e:
		print e
	return True

def SplitInheritencePattern(InpDir, VCF, PED):
	vcfs = utils.get_files(InpDir,'.tsv')
	peds = utils.get_files(InpDir,'.ped')
	vcf_peds = MatchVcfPed(vcfs,peds)
	#vcf_peds = [[VCF,PED]]	
	for vcf,ped in vcf_peds:
		Ped = Pedigree(ped)
		print "Processing vcf: %s\tpedigree: %s" % (vcf,ped)
		fout_AR = open(utils.GetBaseName(vcf)+'_AR.tsv','wb')
		fout_AR.write('Autsomal Recessive Variants\n')
		fout_AD = open(utils.GetBaseName(vcf)+'_AD.tsv','wb')
		fout_AD.write('Autsomal Dominant Variants\n')
		fout_XL = open(utils.GetBaseName(vcf)+'_XL.tsv','wb')
		fout_XL.write('X-linked Variants\n')
		fout_CH = open(utils.GetBaseName(vcf)+'_CH.tsv','wb')
		fout_CH.write('Compound Heterozygote\n')

		fin = open(vcf,'rb')
		header = fin.readline()
		fout_AR.write(header)
		fout_AD.write(header)
		fout_XL.write(header)
		fout_CH.write(header)
		headerList = header.strip().split('\t')
		RecordLen = len(headerList)
		Gene_idx = headerList.index('Annotation')
		Model_idx = headerList.index('GeneticModels')
		Format_idx = headerList.index('FORMAT')	
		CH_buffer = []
		LastGene = None
		for l in fin:
			llist = l.strip().split('\t')
			Gene = llist[Gene_idx]
			Model = llist[Model_idx].split(':')[-1]
			if 'AR' in Model and selectAR(llist, headerList, Ped, Model, Format_idx): 
				fout_AR.write(l)
			if 'AD' in Model and selectAD_DN(llist,headerList, Ped, Model, Format_idx):
				fout_AD.write(l)
			if 'X' in Model:
				fout_XL.write(l)

			if Gene != LastGene:
				if LastGene != None:
					LookCompoundHeter(CH_buffer,fout_CH,RecordLen,Format_idx, headerList, Ped)
				CH_buffer = [l]
				LastGene = Gene
			else:
				CH_buffer.append(l)
		LookCompoundHeter(CH_buffer,fout_CH,RecordLen,Format_idx, headerList, Ped)
		fout_AR.close()
		fout_AD.close()
		fout_XL.close()
		fout_CH.close()
		new_dir = os.path.join(InpDir,utils.GetBaseName(vcf))
		print new_dir
		if not os.path.exists(new_dir):
			os.mkdir(new_dir)
		shutil.copy(vcf, new_dir)
		shutil.copy(ped, new_dir)
		for dest in [utils.GetBaseName(vcf)+'_AR.tsv',utils.GetBaseName(vcf)+'_AD.tsv',utils.GetBaseName(vcf)+'_XL.tsv',utils.GetBaseName(vcf)+'_CH.tsv']:
			if os.path.exists(os.path.join(new_dir, dest)):
				os.remove(os.path.join(new_dir, dest))
			shutil.move(dest, new_dir)

def main():
	DIR, VCF,PED = GetOptions()
	SplitInheritencePattern(DIR, VCF, PED)
	return

if __name__=='__main__':
	main()
