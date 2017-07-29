import sys

def doc_info(sample_id, doc_file, flag_file):
   info = {}
   with open(flag_file) as f:
       for line in f:
           if 'total' in line:
               if 'rawreads' not in info:
                   info['rawreads'] = int(line.split()[0])
           if 'mapped' in line:
               if 'mappedreads' not in info:
                   info['mappedreads'] = int(line.split()[0])
   with open(doc_file) as f:
       head = f.readline().strip().split()
       doc_info = dict(zip(head, f.readline().strip().split()))
       info['targetedbases'] = int(doc_info['total'])
       info['mean'] = float(doc_info['mean'])
       info['median'] = float(doc_info['granular_median'])
       info['DP1'] = float(doc_info['%_bases_above_1'])
       info['DP10'] = float(doc_info['%_bases_above_10'])
       info['DP15'] = float(doc_info['%_bases_above_15'])
   
   info['Mapped_perc'] = 100.0 * info['mappedreads'] / info['rawreads']
   info['Reads_mapped_on_target'] = float(info['targetedbases'])  /  info['mappedreads']
   write_info = [sample_id, info['rawreads'], info['mappedreads'], info['targetedbases'], info['mean'], info['median'],
   info['DP1'], info['DP10'], info['DP15'], info['Mapped_perc'], info['Reads_mapped_on_target'] ]
   return write_info

def get_fname(path):
	return path.split('/')[-1]
def get_sample(fname):
	return fname.split('.')[0]

def make_dict(fin2):
	List = [ line.strip() for line in fin2.readlines()]
	res = {}
	for line in List:
		name = get_sample(get_fname(line))
		res[name] = line.strip()
	return res

def match_samples(sample_summary, flagstat):
	fin1 = open(sample_summary,'rb')
	fin2 = open(flagstat,'rb')
	flag_dict = make_dict(fin2)
	for l in fin1:
		sample = get_sample(get_fname(l))
		#print sample
		if sample in flag_dict:
			
			yield sample, l.strip(), flag_dict[sample]

def WriteSum(sample_summary,flagstat):
	fw = open('result.txt', 'w')
	write_info = ['sample_id', '#Raw_Reads', '#Mapped_Reads', '#Bases_mapped_on_Target', 'Avg_depth', 'Median_depth', 'D1(%)',  'D10(%)', 'D15(%)', '%Mapped', '%Reads_mapped_on_target']
	fw.write('\t'.join(write_info) + '\n')

	#sample_id = '1-14798-01'
	#doc_file, flag_file = '1-14798-01-004-004_S52_L004_001.bwamem.mkdup.DoC.sample_summary', '1-14798-01-004-004_S52_L004_001.bwamem.flagstat'
	for sample_id, doc_file, flag_file in match_samples(sample_summary, flagstat):
		try:	
			write_info = doc_info(sample_id, doc_file, flag_file)
			fw.write('\t'.join(map(str, write_info)) + '\n')
		except ZeroDivisionError:
			print "Error with", sample_id, doc_file, flag_file
	fw.close()

def main():
	sample_summary = sys.argv[1]
	flagstat = sys.argv[2]
	WriteSum(sample_summary, flagstat)

if __name__=='__main__':
	main()
