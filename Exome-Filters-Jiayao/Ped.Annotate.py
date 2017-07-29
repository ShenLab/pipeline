#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Ped.Annotate.py
# Annotation Pedigree file with Relatedness2 and sexcheck file
#=========================================================================

import argparse
import pprint

pp = pprint.PrettyPrinter(indent=4)

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--ped', type=str, help='Pedigree file .ped to be annoteted')
	parser.add_argument('-r', '--relate', type=str, help='Relatedness file .relatedness2 to infer kinship from')
	parser.add_argument('-s', '--sex', type=str, help='Sexcheck file .sexcheck to infer sex from')
	parser.add_argument('-o', '--out', type=str, help='Name of annotated Pedigree file')
	args = parser.parse_args()
	if args.out == None:
		args.out = args.ped.rstrip('ped')+'Annotated.ped'
	return args.ped, args.relate, args.sex, args.out

def GetRelate(RelFil):
	if RelFil != None:
		fin = open(RelFil, 'rb')
		fin.readline()
		res = {}
		for l in fin:
			ind1, ind2, a,b,c,d, phi = l.strip().split('\t')
			if ind1 not in res:
				res[ind1] = {}
			if ind2 not in res[ind1]:
				res[ind1][ind2] = phi
		#pp.pprint(res)
	else:
		res = None
	return res

def GetSex(SexFil):
	res = None
	if SexFil != None:
		fin = open(SexFil, 'rb')
		fin.readline()
		res = {}
		for l in fin:
			#print l
			fid, a, b, snpsex, d, e = l.strip().split()
			#print fid, a, b, snpsex, d, e 
			res[fid] = [snpsex,e]
		#pp.pprint(res)
	return res

def Annotate_2(PedFil, RelateDict, SexDict, OutFil):
	fin = open(PedFil, 'rb')
	fout = open(OutFil, 'wb')
	for l in fin:
		if l.startswith('#'):
			fout.write(l.strip()+'\tRelateness\tSexCheck\n')
			continue
		fam, _id, fa, mo, sex, pheno = l.strip().split('\t')[:6]
		#print fam, _id, fa, mo, sex, pheno
		if fam == _id: #proband
			re_fa = RelateDict[_id][fa]
			re_mo = RelateDict[_id][mo]
			sex = SexDict[_id]
			#fout.write(l.strip()+'\t{}\t{}\n'.format(re_fa+'/'+re_mo, ':'.join(sex)))
			fout.write(l.strip()+'\t{}\t{}\n'.format(re_fa+'/'+re_mo, ':'.join(sex)))
		else:
			re = RelateDict[fam][_id]
			sex = SexDict[_id]
			fout.write(l.strip()+'\t{}\t{}\n'.format(re, ':'.join(sex)))
	return

def Annotate(PedFil, RelateDict, SexDict, OutFil):
	fin = open(PedFil, 'rb')
	fout = open(OutFil, 'wb')
	for l in fin:
		if l.startswith('#'):
			fout.write(l.strip()+'\tRelateness\tSexCheck\n')
			continue
		fam, _id, fa, mo, sex, pheno = l.strip().split('\t')[:6]
		#print fam, _id, fa, mo, sex, pheno
		if fa != '0':
			re_fa = RelateDict[_id][fa]
		else:
			re_fa = '.'
		if mo != '0':
			re_mo = RelateDict[_id][mo]
		else:
			re_mo = '.'
			sex = SexDict[_id]
		fout.write(l.strip()+'\t{}\t{}\n'.format(re_fa+'/'+re_mo, ':'.join(sex)))
	return

def main():
	PedFil, RelFil, SexFil, OutFil = GetOptions()
	RelateDict = GetRelate(RelFil)
	SexDict = GetSex(SexFil)
	Annotate(PedFil, RelateDict, SexDict, OutFil)
	return


if __name__ == '__main__':
	main()
