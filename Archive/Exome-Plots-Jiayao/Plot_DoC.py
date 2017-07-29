#!/home/local/users/jw/bin/python2.7
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
#
#=========================================================================

from optparse import OptionParser
from utils import *
import matplotlib as mpl
mpl.use('Agg')
#plt.use('qt5agg')
import matplotlib.pyplot as plt


def GetOptions():
	parser = OptionParser()
	parser.add_option('-d', '--dir', dest='DIR', metavar='DIR',
			help='Dirctory contains sample summary')
	parser.add_option('-m', '--mode', dest='MODE', metavar='MODE',
			help='Mode to excute, 1, prepare input table; 2, plot D10-D15/Mean plot')
	parser.add_option('-t', '--table', dest='Table',
			metavar='Table', help='DP Table to draw plot from')
	parser.add_option('-o', '--out', dest='Outname',
			metavar='Outname', help='Outname')
	(options, args) = parser.parse_args()

	return options.MODE, options.DIR, options.Table, options.Outname


def Draw_Table(DIR, outname):
	sample_file = get_files(DIR, 'DoC.sample_summary')
	fout = open(outname, 'wb')
	headers = ['sample_id', 'total', 'mean', 'granular_third_quartile', 'granular_median', 'granular_first_quartile',
			'%_bases_above_1', '%_bases_above_5', '%_bases_above_10', '%_bases_above_15', '%_bases_above_20']
	fout.write('\t'.join(headers) + '\n')
	for sample in sample_file:
		with open(sample) as fin:
			data = fin.readlines()[1]
			fout.write(data)
	fout.close()


def Plot_D15(Table):
	fin = open(Table, 'rb')
	fin.readline()
	Means, D10, D15 = [], [], []
	for l in fin:
		tmps = l.strip().split('\t')
		if float(tmps[8]) > 80:
			Means.append(float(tmps[2]))
			D10.append(float(tmps[8]))
			D15.append(float(tmps[9]))
	fin.close()
	Plot_scatter(Means, D10, D15)


def Plot_scatter(Means, D10, D15):
	# print Means
	# print Means,D10,D15
	# D10:
	# print "Plot"
	plt.scatter(Means, D10, color='b', marker='o', label='D10', alpha=1, s=2)
	plt.scatter(Means, D15, color='r', marker='o', label='D15', alpha=1, s=2)
	plt.xlabel('Mean Coverage')
	plt.ylabel('% Targeted Bases')
	plt.title('Mean-DP15')
	# plt.axis([0, 201, 0, 0.05])ppppp
	#plt.xticks(np.arange(0, 501, 10),rotation='vertical')
	plt.legend(loc='upper right')
	plt.grid(True)
	#plt.show()
	plt.savefig('DoC.png')


def main():
	mode, DIR, Table, outname = GetOptions()
	if mode == '1':
		Draw_Table(DIR, outname)
	elif mode == '2':
		Plot_D15(Table)
	else:
		if outname == None:
			outname, Table = "DoC.txt", "DoC.txt"
		Draw_Table(DIR, outname)
		Plot_D15(Table)

if __name__ == '__main__':
	main()
