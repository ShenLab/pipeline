#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# SplitBed.py
#========================================================================================================

import argparse
import gzip

class Interval():
    def __init__(self, CHR, START, END):
        self.CHR = CHR
        #self.START=max(0, int(START)-10)
        #self.END=int(END) + 10
        self.START = int(START)
        self.END = int(END)
        self.LENGTH = int(self.END) - int(START) + 1

class SplitBed:
	def __init__(self, targetFil, N):
		self.targetFil = target
		self.N = N
	def run(self):
		# Read All intervals into list
		fin = open(self.targetFil, 'rb')
		intervals = []
		Total_length = 0
		for l in fin:
			if l.startswith('#'):
				continue
			CHR, START, END = l.strip().split()[:3] 
			CHR = CHR.strip('chr') 
			interval = Interval(CHR, START, END)
			intervals.append(interval)
			Total_length += interval.LENGTH
		window_length = ceil((Total_length + 0.0) / self.N)

		current_window_len = 0
		current_window = []
		num_of_bed = 0
		bedfiles = []
		for interval in intervals:
			current_window.append(interval)
			current_window_len += interval.LENGTH 
			if current_window_len >= window_length:
				current_window_len = 0
				bed_name = "Splited_" + str(num_of_bed) + '.bed'
				num_of_bed += 1
				make_bed(bed_name, current_window)
				bedfiles.append(bed_name)
				current_window = []
		bed_name = "New_bed_" + str(num_of_bed) + '.bed'
		if len(current_window) > 0:
			make_bed(bed_name, current_window)
			bedfiles.append(bed_name)

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-t','--target', type=str, help = 'Target Bed file to be split')
	parser.add_argument('-n','--Num', type=str, help = 'How many parts you want to split?')
	args = parser.parse_args()
	
	return args.target, args.Num

def main():
	target, Num = GetOptions()
	ins = SplitBed(target, Num)
	ins.run()
	return

if __name__=='__main__':
	main()
