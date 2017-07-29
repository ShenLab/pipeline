#!/home/yufengshen/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# csv2ped.py
#========================================================================================================

import argparse
import csv
import re

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--csv', type=str, help = 'csv file contains pedigree information')
    parser.add_argument('-p','--ped', type=str, help = 'Pedigree file output the result')
    args = parser.parse_args()
    return args.csv, args.ped

class Fam():
    def __init__(self, FAMID):
        self.FamID = FAMID
        self.Member = []
    def AddProband(self, Proband):
        self.Proband = Proband
    def AddMember(self, Sample):
        self.Member.append(Sample)
    def FormRecord(self):
        res = [self.Proband.FormRecord()]
        for member in self.Member:
            res.append(member.FormRecord())
        return ''.join(res)

class Sample():
    def __init__(self, FAM, ID, comment):
        self.FamID = FAM
        self.ID = ''.join(ID.split())
        self.FatherID = '0'
        self.MotherID = '0'
        self.Gender = '0'
        self.Phenotype = '0'
        self.Comment = comment
    def AddGender(self, Gender):
        self.Gender = Gender
    def AddFather(self, FatherID):
        self.FatherID = FatherID
    def AddMother(self, MotherID):
        self.MotherID = MotherID
    def FormRecord(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.FamID, self.ID, self.FatherID, self.MotherID, self.Gender, self.Phenotype)

def ReadCSV(CSV):
    res = []
    fin = open(CSV, 'rb')
    csv_reader = csv.reader(fin)
    next(csv_reader, None)
    for row in csv_reader:
        res.append(row)
    fin.close()
    return res

def FormPed(csv_list, PED):
    # INIT FAM at first around
    FAMs = {}
    for row in csv_list:
        sample, _id, fam, relationship, phenotype = row
        Individual = Sample(fam, sample, [relationship.strip(), phenotype.strip()])
        if Individual.Comment[1].lower() == 'unaffected':
            Individual.Phenotype = '1'
        else:
            Individual.Phenotype = '2'
        if fam not in FAMs:
            FAMs[fam] = Fam(fam)
        #if Individual.Comment[0].lower() == 'proband':
        #    FAMs[fam].AddProband(Individual)
        #if:
        FAMs[fam].AddMember(Individual)
    # Format the Fam 
    fout = open(PED, 'wb')
    fout.write('#FAMID\tProbandID\tFatherID\tMotherID\tGender\tPhenotype\tComment1\tComment2\n')
    FAM_list = sorted(FAMs.values(), key = lambda x:x.FamID)
    for FAM in FAM_list:
        for Member in FAM.Member:
            if Member.Comment[0].lower() == 'proband':
                Proband = Member
                Proband.FamID = Proband.ID
        for Member in FAM.Member:
            if re.match('father',Member.Comment[0].lower()):
                Father = Member
                Father.FamID = Proband.ID
                Father.AddGender('1')
                Proband.FatherID = Father.ID
        for Member in FAM.Member:
            if re.match('mother',Member.Comment[0].lower()):
                Mother = Member
                Mother.FamID = Proband.ID
                Mother.AddGender('2')
                Proband.MotherID = Mother.ID
        fout.write(Proband.FormRecord())
        fout.write(Father.FormRecord())
        fout.write(Mother.FormRecord())
def main():
    CSV, PED = GetOptions()
    csv_list = ReadCSV(CSV)
    FormPed(csv_list, PED)
    return

if __name__=='__main__':
    main()
