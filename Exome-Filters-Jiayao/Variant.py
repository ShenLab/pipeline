# Contains several class for parse the data need to process vcf files.
import re

# =================================================================================
# The sample in a vcf record.
# =================================================================================


class Sample():
    # 1-00009
    # GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
    def __init__(self, name, Format, data):
        self.name = name  # Sample name 1-00009
        # Format ["GT", "AD", "DP", "GQ", "PL"]
        self.format = Format.split(':')
        # Data ["0/0", "7,0", "7", "18", "0,18,270"]
        self.data = data.split(":")
        self.dict = {}
        for i, item in enumerate(self.format):
            self.dict[item] = self.data[i]
            # If sample don't have genotpye, no need to parse further
            if self.dict['GT'] == './.':
                break
        # GT '0/0' --> ['0', '0']
        self.GT = re.findall('[\d.]', self.dict['GT'])

    def get_format(self):
        return ":".join(self.format)

    def get_data(self):
        return ":".join(self.data)
# =================================================================================

# =================================================================================
# INFO Field of a VCF Record
# Consist of two parts, from caller & from annovar
# For INFO from annovar, there are (Num of alt allele) sets of INFOs
# =================================================================================


class INFO():
    def __init__(self, alts, raw_string):
        # This INFO group is unique for each site. Come from variants caller or
        # other annotations.
        self.get_unique(raw_string)
        # This INFO group is unique for each alt allele, Come from ANNOVAR
        # annotation that give each allele one annotation.
        self.get_annovar(alts, raw_string)

    # A set of commen annotation is before the first "ANNOVAR_DATE" and after
    # last "ALLELE_END"
    def get_unique(self, raw_string):
        self.raw_string = raw_string
        unique = re.compile("(.+)ANNOVAR_DATE")
        groups = unique.findall(raw_string)
        # print "unique info:",groups[0]
        tmp = (groups[0]).split(';')
        self.unique = {}
        for item in tmp:
            kv = item.split('=')
            if len(kv) == 2:
                k, v = kv
                self.unique[k] = v
            elif len(kv) == 1:
                # print 'Without Value:',kv
                self.unique[kv[0]] = ""
            else:
                print "Unexpected info:", kv

    # A set of annovar annotation is between "ANNOVAR_DATE" and "ALLELE_END"
    def get_annovar(self, alts, raw_string):
        self.annovar = []
        get_allele_anno = re.compile("ANNOVAR_DATE=[\d-]+;.+?ALLELE_END")
        groups = get_allele_anno.findall(raw_string)
        # print "annovar info:",groups
        for i, alt in enumerate(alts):
            tmp = {}
            for item in groups[i].split(';'):
                kv = item.split('=')
                if len(kv) == 2:
                    k, v = kv
                    tmp[k] = v
                elif len(kv) == 1:
                    # print 'Without Value:',kv
                    tmp[kv[0]] = ""
                else:
                    print "Unexpected info:", kv
            self.annovar.append(tmp)


# =================================================================================
# A Variant record in vcf file
# record is the line of the variants in vcf file
# headers is the header line (startswith '#'), contains what info and samples this variants should have
# =================================================================================
class Variant():
    def __init__(self, record, headers):
        self.headers = headers
        self.record = record.strip().split('\t')
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = self.record[0:9]
        self.Chrom = CHROM
        self.Pos = POS
        self.Id = ID
        self.Ref = REF
        self.Alt = ALT.split(',')  # A list
        self.Qual = QUAL
        self.Filter = FILTER
        self.Info_str = INFO
        self.GetInfo(INFO)
        self.Format = FORMAT
        self.GetSample(record)
    # Map the informations into a dictionary.

    def GetInfo(self, INFO_string):
        self.Info = INFO(self.Alt, INFO_string)
    # Generate sample class for each sample in the variants record.

    def GetSample(self, record):
        self.Sample = {}
        for i, sample in enumerate(self.headers[9:]):
            # print headers
            # print record
            self.Sample[sample] = Sample(
                sample, self.Format, self.record[i + 9])
    # Form the Info str from a given field sequence.

    def write_info(self, info_seq=None):
        self.Info_str = self.Info.raw_string
        return self.Info.raw_string

    # form back vcf record with first 8 fileds.
    def write_basic(self):
        return '\t'.join([self.Chrom, self.Pos, self.Id, self.Ref, ','.join(self.Alt), self.Qual, self.Filter, self.Info_str, self.Format])
    # form back vcf record with samples given.

    def write_samples(self, samples):
        res = []
        try:
            for i, sample in enumerate(samples):
                res.append(':'.join(self.Sample[sample].data))
            # print res
            return '\t'.join(res)
        except KeyError:
            print "No sample %s in this VCF file !" % sample
            exit()
    # Combine basic and samples out put the final result line.

    def write(self, samples):
        return self.write_basic() + '\t' + self.write_samples(samples) + '\n'

    # print samples
    def show_samples(self):
        for k, v in self.Sample.items():
            print '%s:%v' % (k, v)

    # Fet a unique variant key. key is chr-pos-ref-alt
    def get_key(self):
        return '|'.join([self.Chrom, self.Pos, self.Ref, ','.join(self.Alt)])

    # Get a unique variant-sample key. key is chr-pos-allele1-allele2.
    def get_sample_key(self, sample_name):
        alleles = [self.Ref]
        alleles.extend(self.Alt)
        GT = sorted(self.Sample[sample_name].GT)
        allele1 = alleles[int(GT[0])]
        allele2 = alleles[int(GT[1])]
        return '|'.join([self.Chrom, self.Pos, allele1, allele2])

# =================================================================================
# Parse a fam in Ped file.
# fam    proband    father    mother
# out_name proband_father_mother.vcf
# =================================================================================


class Pedigree():
    def __init__(self, line):
        self.fam, self.proband, self.father, self.mother = line.strip().split('\t')[
            0:4]
        self.out_name = '_'.join(
            [self.proband, self.father, self.mother]) + '.vcf'


def get_pedigrees(ped_file):
    fin = open(ped_file, 'rb')
    res = {}
    fin.readline()
    for line in fin:
        pedigree = Pedigree(line)
        if pedigree.fam in res or pedigree.proband == '0':
            continue
        else:
            res[pedigree.fam] = pedigree
    return res.values()
