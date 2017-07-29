#!/home/local/users/jw/bin/python2.7
# =================================================================================
# Script used for Seperate a vcf file into n vcf files according to a pedifree file
# =================================================================================
#from Filters import *
from Variant import *
import argparse
import sys
import os
import re
import gzip

# =================================================================================
# Get the basename of a vcf file, should be name of a variant caller
# E.g: fatk, st, pt, fb
def get_basename(name):
	return name.split('/')[-1].split('.')[0]

# =================================================================================
# Trim header line with a trios.
# the 8 basics + proband father mother
# =================================================================================
def trim_head(header,fam):
	if fam.father != '0' and fam.mother != '0':
		return '\t'.join(header[0:9])+'\t'+'\t'.join([fam.proband,fam.father,fam.mother])+'\n'
	elif fam.father != '0':
		return '\t'.join(header[0:9])+'\t'+'\t'.join([fam.proband,fam.father]) + '\n'
	elif fam.mother != '0':
		return '\t'.join(header[0:9])+'\t'+'\t'.join([fam.proband,fam.mother]) + '\n'
	else:
		return '\t'.join(header[0:9])+'\t'+'\t'.join([fam.proband]) + '\n'
		
# =================================================================================
# Seperate a vcf file by a ped_file
# Given a work_dir, fout splited files here and put into a new dir with associate name.
# =================================================================================
def seperate_by_trios(work_dir,vcf_file,fams,filters=None):
	if not work_dir.endswith('/'):
		work_dir += '/'
	os.chdir(work_dir)
	vcf_hand = gzip.open(vcf_file,'rb')
	vcf_file_basename = get_basename(vcf_file)
	print "Processiing",vcf_file_basename,"variants"
	#open len(fams) handle to save variants in each proband
	sub_vcf_names = []
	sub_vcf_hands = []
	for fam in fams:
		sub_name = vcf_file_basename+'_'+fam.out_name
		sub_vcf_names.append(sub_name)
		sub_vcf_hands.append(open(sub_name,'wb'))
	meta_info = []
	for line in vcf_hand:
		#Skip headers
		if line.startswith('##'):
			meta_info.append(line)
		elif line.startswith('#'):
			header = line.strip().split('\t')
			for fam_i,fam in enumerate(fams):
				sub_vcf_hands[fam_i].write(''.join(meta_info))
				sub_vcf_hands[fam_i].write(trim_head(header,fam))
		#read variants
		else: 
			var = Variant(line,header)
			for fam_i,fam in enumerate(fams):
				# If on this site, proban have an un-wild genotype and parients have a clear genotype, write file
				#print fam.proband
				"""
				#works for compete trios, not for incomplete trios.
				try:
					if (var.Sample[fam.proband].GT != ['0','0']) and ('.' not in var.Sample[fam.proband].GT) and ('.' not in var.Sample[fam.father].GT) and ('.' not in var.Sample[fam.mother].GT):
						#print var.Sample[fam.proband].GT
						if filters == None:	
							sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.father,fam.mother]))
						else:
							if pass_filter(var,fam,vcf_file_basename):
								sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.father,fam.mother]))
				except KeyError as e:
					print e
					print "No sample:" ,fam.proband
					pass
				"""
				# Try to make it work for incomplete trios
				try:
					#complte trios
					if fam.proband in var.Sample and fam.father in var.Sample and fam.mother in var.Sample:
						if (var.Sample[fam.proband].GT != ['0','0']) and ('.' not in var.Sample[fam.proband].GT) and ('.' not in var.Sample[fam.father].GT) and ('.' not in var.Sample[fam.mother].GT):
							if filters == None:	
								sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.father,fam.mother]))
							else:
								if pass_filter(var,fam,vcf_file_basename):
									sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.father,fam.mother]))
					elif fam.proband in var.Sample and fam.father in var.Sample:
						if (var.Sample[fam.proband].GT != ['0','0']) and ('.' not in var.Sample[fam.proband].GT) and ('.' not in var.Sample[fam.father].GT):
							if filters == None:	
								sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.father]))
							else:
								if pass_filter(var,fam,vcf_file_basename):
									sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.father]))
						
					elif fam.proband in var.Sample and fam.mother in var.Sample:
						if (var.Sample[fam.proband].GT != ['0','0']) and ('.' not in var.Sample[fam.proband].GT) and ('.' not in var.Sample[fam.mother].GT):
							if filters == None:	
								sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.mother]))
							else:
								if pass_filter(var,fam,vcf_file_basename):
									sub_vcf_hands[fam_i].write(var.write([fam.proband,fam.mother]))

					elif fam.proband in var.Sample:
						if (var.Sample[fam.proband].GT != ['0','0'] and ('.' not in var.Sample[fam.proband].GT)): 
							if filters == None:	
								sub_vcf_hands[fam_i].write(var.write([fam.proband]))
							else:
								if pass_filter(var,fam,vcf_file_basename):
									sub_vcf_hands[fam_i].write(var.write([fam.proband]))

					else:
						print fam.proband,"Not in vcf file, please check the individual ids"
						exit()
				except KeyError as e:
					print e
					print line
					exit()
					pass
	vcf_hand.close()
	for sub_hand in sub_vcf_hands:
		sub_hand.close()

	# mv all new vcf to right place
	current_sub_dir = work_dir+vcf_file_basename+'_splited_trios'
	if not os.path.exists(current_sub_dir):
		os.makedirs(current_sub_dir)
	for sub_name in sub_vcf_names:
		os.rename(sub_name, current_sub_dir+'/'+sub_name)
		#print current_sub_dir+'/'+sub_name
	print "Done!"

# =================================================================================
# Given a ped file, a list of vcf and a work dir.
# Seperate each vcf file according to a pedigree file. output splited vcf files 
# =================================================================================
def Seperate(work_dir,VCFs,ped_file,debug):
	# Example:
	#	work_dir = "/home/local/users/jw/Consensus_Pipeline/trios_50/Analysis/SNV/"
	#	gatk_snv = "gatk.snv.vcf"
	#	st_snv = "st.snv.vcf"
	#	pt_snv = "pt.snv.vcf"
	#	ped_file = "/home/local/users/jw/Consensus_Pipeline/trios_50/trios_50.ped"
	work_dir = os.path.abspath(work_dir)
	print "Work Dir is:",work_dir
	fams = get_pedigrees(ped_file)
	for fam in fams:
		print fam.proband,fam.father,fam.mother
	for vcf in VCFs:
		seperate_by_trios(work_dir,vcf,fams,filters=['GT','ExAC_All'])
	return
# =================================================================================

# =================================================================================
# Given a variant and trio
# Filter the variant for some contitions, like GQ, ExAC_All
# =================================================================================
def pass_filter(var,fam,vcf_file_basename):
	if vcf_file_basename == 'gatk':
		if float(var.Sample[fam.proband].dict['GQ']) <= 20:
			return False
	elif vcf_file_basename == 'st':
		if pick_PL(var.Sample[fam.proband].dict['PL']) <= 20:
			return False
	elif vcf_file_basename == 'pt':
		if float(var.Sample[fam.proband].dict['GQ']) <= 20:
			return False
	elif vcf_file_basename == 'fb':
		if pick_GL(var.Sample[fam.proband].dict['GL']) >= -10:
			return False
	else:
		try:
			if float(var.Sample[fam.proband].dict['GQ']) <= 20:
				return False
		except KeyError:
			print var.Sample[fam.proband].dict
			exit()
	try:
		if float(var.Info.annovar[int(var.Sample[fam.proband].GT[1])-1]['ExAC_ALL']) > 0.01:
			return False
	except ValueError:
		return True
	return True
def pick_PL(PL):
	PL = map(float,PL.split(','))
	PL.sort()
	return PL[1]
def pick_GL(GL):
	GL = map(float,GL.split(','))
	GL.sort()
	return GL[-2]
# =================================================================================
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--dir", type=str,required=True,
        help="<Required> Enter the work_dir to generate the results")
	parser.add_argument("-v","--vcf_file",type=str,nargs='+',required=True,
		help="<Required> Enter the vcf file you wanna seperate. E.g: gatk.vcf")
	parser.add_argument('-p',"--pedigree",type=str,default='/home/local/users/jw/Consensus_Pipeline/trios_50/trios_50.ped',
		help="<Required> pedigree file you want use for seperated by")
	parser.add_argument('--debug',type=int,default=0,choices=[0,1],
		help="Turn on Debug mod, output many mid results for debugging default=0")
	args = parser.parse_args()
	print args.vcf_file
	Seperate(args.dir,args.vcf_file,args.pedigree,args.debug)
# =================================================================================
if __name__=='__main__':
	main()
