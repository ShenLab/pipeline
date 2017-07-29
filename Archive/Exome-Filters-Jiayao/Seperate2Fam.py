#!/home/local/users/jw/bin/python2.7
# =================================================================================
# Script used for Seperate a vcf file into n vcf files according to a pedifree file
# =================================================================================
from Filters import *
from Variant import Variant
import argparse
import sys
import time
import os
import re
import gzip
# =================================================================================
# Manipulate Pedgrees


def WriteSeperatePedFile(fam, PedHeader, buff):
    fout = open(fam + '.ped', 'wb')
    fout.write(PedHeader)
    for l in buff:
        fout.write(l)
    fout.close()


def ReadPedigreeFromPedFile(PedFile):
    fin = open(PedFile)
    PedHeader = fin.readline()
    Res = []  # The list store pedigrees
    buff = []  # Used to store lines from ped file in same family
    Last_fam = None
    counter = 0
    for l in fin:
        llist = l.strip().split('\t')
        fam = llist[0]
        if fam != Last_fam:
            if Last_fam != None:
                Res.append(Pedigree(buff))
                WriteSeperatePedFile(Last_fam, PedHeader, buff)
                counter += 1
            buff = [l]
            Last_fam = fam
        else:
            buff.append(l)
    Res.append(Pedigree(buff))
    WriteSeperatePedFile(Last_fam, PedHeader, buff)
    counter += 1
    print "Total %d Trio find" % counter
    return Res


class Pedigree():
    def __init__(self, buff):
        self.famid = buff[0].split('\t')[0]
        self.out_name = self.famid + '.vcf'
        self.individuals = []
        self.ids = []
        for l in buff:
            llist = l.strip().split('\t')
            # (SampleID, Phenotype)
            self.individuals.append((llist[1], llist[5]))
            self.ids.append(llist[1])

    def show(self):
        print "\nFamily:", self.famid
        for k, v in self.individuals:
            print "SampleID: %s\tPhenotype: %s" % (k, v)
        print

# =================================================================================

# =================================================================================
# Get the basename of a vcf file, should be name of a variant caller
# E.g: fatk, st, pt, fb


def get_basename(name):
    return name.split('/')[-1].split('.')[0]

# =================================================================================
# Trim header line with a trios.
# the 8 basics + proband father mother
# =================================================================================


def trim_head(header, indiList):
    return '\t'.join(header[0:9]) + '\t' + '\t'.join(indiList) + '\n'

# =================================================================================
# Seperate a vcf file by a ped_file
# Given a work_dir, fout splited files here and put into a new dir with associate name.
# =================================================================================


def seperate_by_fam(work_dir, vcf_file, fams, filters=None):
    if not work_dir.endswith('/'):
        work_dir += '/'
    os.chdir(work_dir)
    if vcf_file.endswith('.vcf.gz'):
        vcf_hand = gzip.open(vcf_file, 'rb')
    else:
        vcf_hand = open(vcf_file, 'rb')
    vcf_file_basename = get_basename(vcf_file)
    print "Processiing", vcf_file_basename, "variants"
    # open len(fams) handle to save variants in each proband
    sub_vcf_names = []
    sub_vcf_hands = []
    for fam in fams:
        #sub_name = vcf_file_basename + '_' + fam.out_name
        sub_name = fam.out_name
        sub_vcf_names.append(sub_name)
        sub_vcf_hands.append(open(sub_name, 'wb'))
    meta_info = []
    counter = 1
    s_time = time.time()
    for line in vcf_hand:
        # Skip headers
        if line.startswith('##'):
            meta_info.append(line)
        elif line.startswith('#'):
            header = line.strip().split('\t')
            for fam_i, fam in enumerate(fams):
                sub_vcf_hands[fam_i].write(''.join(meta_info))
                sub_vcf_hands[fam_i].write(trim_head(header, fam.ids))
        # read variants
        else:
            if counter %1000 == 0:
                print "Read %d records, in %.3fs" %(counter, time.time() - s_time)
                s_time = time.time()
            var = Variant(line, header)
            for fam_i, fam in enumerate(fams):
                # If on this site, proband have an un-wild genotype and parients have a clear genotype, write file
                # print fam.proband
                flag = False
                for k, v in fam.individuals:
                    if v == '2':  # Phenotype
                        if var.Sample[k].GT != ['0', '0'] and ('.' not in var.Sample[k].GT):
                            flag = True
                if flag:
                    if filters == None:
                        sub_vcf_hands[fam_i].write(var.write(fam.ids))
                    else:
                        if pass_filter(var, fam, vcf_file_basename):
                            sub_vcf_hands[fam_i].write(var.write(fam.ids))
        counter += 1
    vcf_hand.close()
    for sub_hand in sub_vcf_hands:
        sub_hand.close()

    # mv all new vcf to right place
    current_sub_dir = work_dir + vcf_file_basename + '_splited_trios'
    if not os.path.exists(current_sub_dir):
        os.makedirs(current_sub_dir)
    for sub_name in sub_vcf_names:
        ped_name = sub_name.rstrip('vcf')+'ped'
        print sub_name, sub_name
        os.rename(sub_name, current_sub_dir + '/' + sub_name)
        os.rename(ped_name, current_sub_dir + '/' + ped_name)
        # print current_sub_dir+'/'+sub_name
    print "Done!"

# =================================================================================
# Given a ped file, a list of vcf and a work dir.
# Seperate each vcf file according to a pedigree file. output splited vcf files
# =================================================================================


def Seperate(work_dir, VCFs, ped_file, debug):
    work_dir = os.path.abspath(work_dir)
    fams = ReadPedigreeFromPedFile(ped_file)
    for i,fam in enumerate(fams):
        print i+1
        fam.show()
    for vcf in VCFs:
        seperate_by_fam(work_dir, vcf, fams, filters=['GT', 'ExAC_All'])
        pass
    return
# =================================================================================

# =================================================================================
# Given a variant and trio
# Filter the variant for some contitions, like GQ, ExAC_All
# =================================================================================


def pass_filter(var, fam, vcf_file_basename):
    try:
        flag = False
        for k, v in fam.individuals:
            if v == '2':
                if float(var.Info.annovar[int(var.Sample[k].GT[1]) - 1]['ExAC_ALL']) > 0.01:
                    flag = flag or False
                else:
                    flag = flag or True
        if not flag:
            return False
    except ValueError:
        return True
    return True


def pick_PL(PL):
    PL = map(float, PL.split(','))
    PL.sort()
    return PL[1]


def pick_GL(GL):
    GL = map(float, GL.split(','))
    GL.sort()
    return GL[-2]
# =================================================================================


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", type=str,
                        help="Enter the work_dir to generate the results")
    parser.add_argument("-v", "--vcf_file", type=str, nargs='+', required=True,
                        help="<Required> Enter the vcf file you wanna seperate. E.g: gatk.vcf")
    parser.add_argument('-p', "--pedigree", type=str, default='/home/local/users/jw/Consensus_Pipeline/trios_50/trios_50.ped',
                        help="<Required> pedigree file you want use for seperated by")
    parser.add_argument('--debug', type=int, default=0, choices=[0, 1],
                        help="Turn on Debug mod, output many mid results for debugging default=0")
    args = parser.parse_args()
    if args.dir == None:
        args.dir = os.getcwd()
    Seperate(args.dir, args.vcf_file, args.pedigree, args.debug)


# =================================================================================
if __name__ == '__main__':
    main()
