#!/home/local/users/jw/bin/python2.7
from random import shuffle
from math import ceil


def read_variables():
    f_var = open("VC_freebayes.sh")
    n_cpu = 1
    freebayes, ref_genome, bedfile, bamlist, bcftools = [None] * 5
    for l in f_var:
        if l.startswith('#'):
            continue
        if l.split('=')[0].strip() == 'FREEBAYES':
            freebayes = l.strip().split('=')[1]
        if l.split('=')[0].strip() == 'REF_GENOME':
            ref_genome = l.strip().split('=')[1]
        if l.split('=')[0].strip() == 'TARGET_INTERVAL':
            bedfile = l.strip().split('=')[1]
        if l.split('=')[0].strip() == 'BAM_LIST':
            bamlist = l.strip().split('=')[1]
        if l.split('=')[0].strip() == 'nCPU':
            n_cpu = l.strip().split('=')[1]
    # check required variable
    check_flag = True
    if freebayes == None:
        check_flag = False
        print "Freebayes path not provided."
    if ref_genome == None:
        check_flag = False
        print "Reference Genome Path not provided"
    if bamlist == None:
        check_flag = False
        print "bamlist file not provided"
    if bedfile == None:
        check_flag = False
        print "Target Interval bed file not provided"
    if check_flag == False:
        print "Please provided required files"
        exit()
    print "N0.CPU: %s\nFreebayes path: %s\nBcftools path: %s\nBam list file: %s\nBedfile: %s\nRef genome: %s" % (n_cpu, freebayes, bcftools, bamlist, bedfile, ref_genome)
    return int(n_cpu), freebayes, bamlist, bedfile, ref_genome


# Create new bed files according to how many process to run parallelly and
# target bedfile.
def calculate_bed(n_cpu, bedfile):
    class Interval():
        def __init__(self, CHR, START, END):
            self.CHR = CHR
            self.START = START
            self.END = END
            self.LENGTH = int(END) - int(START) + 1
    f_bed = open(bedfile, 'rb')
    intervals = []
    Total_length = 0
    for l in f_bed:
        CHR, START, END = l.strip().split()[:3]
        CHR = CHR.strip('chr')
        interval = Interval(CHR, START, END)
        intervals.append(interval)
        Total_length += interval.LENGTH
    shuffle(intervals)
    window_length = ceil((Total_length + 0.0) / n_cpu)
    print Total_length, n_cpu, window_length

    current_window_len = 0
    current_window = []
    num_of_bed = 0
    bedfiles = []
    for interval in intervals:
        current_window.append(interval)
        current_window_len += interval.LENGTH
        if current_window_len >= window_length:
            current_window_len = 0
            bed_name = "New_bed_" + str(num_of_bed) + '.bed'
            num_of_bed += 1
            make_bed(bed_name, current_window)
            bedfiles.append(bed_name)
            current_window = []
    bed_name = "New_bed_" + str(num_of_bed) + '.bed'
    if len(current_window) > 0:
        make_bed(bed_name, current_window)
        bedfiles.append(bed_name)
    return bedfiles


def add_zero(CHR):
    if CHR in ['1', '2', '3', '4', '5', '6', '7', '8', '9']:
        return '0' + CHR
    else:
        return CHR


def make_bed(bed_name, current_window):
    current_window.sort(key=lambda x: (add_zero(x.CHR), int(x.START)))
    fout = open(bed_name, 'wb')
    for item in current_window:
        fout.write(item.CHR + '\t' + item.START + '\t' + item.END + '\n')
    fout.close()


def make_script(bedfiles, freebayes, bamlist, ref_genome):
    print "bedfiles:", bedfiles
    script = open('run_freebayes.sh', 'wb')
    script.write('#!/bin/bash\n')
    script.write('FREEBAYES=' + freebayes + '\n')
    script.write('BAMLIST=' + bamlist + '\n')
    script.write('REFGENOME=' + ref_genome + '\n\n')
    script.write('#Run freebayes call variant for each bed file.\n')
    for i, bed in enumerate(bedfiles):
        script.write('nohup $FREEBAYES -f $REFGENOME -t ' + bed +
                     ' -L $BAMLIST 2>stderr_' + str(i) + '.txt > freebayes.' + str(i) + '.vcf &\n\n')


def main():

    n_cpu, freebayes, bamlist, bedfile, ref_genome = read_variables()
    new_beds = calculate_bed(n_cpu, bedfile)
    make_script(new_beds, freebayes, bamlist, ref_genome)


if __name__ == '__main__':
    main()
