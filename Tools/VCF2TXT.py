#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# VCF2TXT.py
#========================================================================================================

import argparse
import csv


def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v','--vcf', type=str, help = 'InputVCF')
	parser.add_argument('-o','--out', type=str, help = 'OutPut Table')
	args = parser.parse_args()
	return args.vcf, args.out

def FlatenVCF(vcf, out):


def main():
	vcf, out = GetOptions()
	FlatenVCF(vcf, out)
	return

if __name__=='__main__':
	main()
