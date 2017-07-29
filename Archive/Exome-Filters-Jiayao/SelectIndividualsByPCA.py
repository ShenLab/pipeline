#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# SelectIndividualsByPCA.py
#========================================================================================================

import argparse

class Individual():
	def __init__(self, llist):
		self.SampleID = llist[0].strip().split(':')[0].strip()
		self.PC1 = float(llist[1].strip())
		self.PC2 = float(llist[2].strip())
		self.PC3 = float(llist[3].strip())
	def Filter(self, egv1, egv2, egv3):
		FlagPC1 = self.PC1 > egv1[0] and self.PC1 < egv1[1]
		FlagPC2 = self.PC2 > egv2[0] and self.PC2 < egv2[1]
		FlagPC3 = self.PC3 > egv3[0] and self.PC3 < egv3[1]
		return FlagPC1 and FlagPC2 and FlagPC3

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v','--vcf', type=str, help = 'VCF file to select individuals from')
	parser.add_argument('-p','--pca', type=str, help = 'PCA result, should have .evec suffix')
	parser.add_argument('-egv1','--egv1', type=str, help = 'Cutoff on eigval1, should be \'low,upper\'')
	parser.add_argument('-egv2','--egv2', type=str, help = 'Cutoff on eigval2, should be \'low,upper\'')
	parser.add_argument('-egv3','--egv3', type=str, help = 'Cutoff on eigval3, should be \'low,upper\'')
	parser.add_argument('-o','--outname', type=str, help = 'Output File Name')
	args = parser.parse_args()
	egv1 = map(float,args.egv1.split(','))
	egv2 = map(float,args.egv2.split(','))
	egv3 = map(float,args.egv3.split(','))
	return args.vcf, args.pca, egv1, egv2, egv3, args.outname

def GetInividualListFromPCAresult(pca, egv1, egv2, egv3, OutName):
	PCA_hand = open(pca, 'rb')
	if OutName == None:
		OutName = pca.split('/')[-1] + '.selected.list'
	Out_hand = open(OutName, 'wb')
	for l in PCA_hand:
		llist = [ x.strip() for x in l.strip().split('    ')]
		if llist[0].startswith('#'):
			continue
		else:
			Indi = Individual(llist)
			if Indi.Filter(egv1, egv2, egv3):
				Out_hand.write(Indi.SampleID + '\n')
	PCA_hand.close()
	Out_hand.close()


def main():
	vcf, pca, egv1, egv2, egv3, outname= GetOptions()
	GetInividualListFromPCAresult(pca, egv1, egv2, egv3, outname)
	return

if __name__=='__main__':
	main()
