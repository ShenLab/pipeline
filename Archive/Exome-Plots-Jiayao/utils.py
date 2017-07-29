#!/home/local/users/jw/bin/python2.7
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# Functions that used many times
#========================================================================================================
import os
from optparse import OptionParser

def get_files(dirc,extension):
	files = []
	print dirc
	for f in os.listdir(dirc):
		if f.endswith(extension):
			#print os.path.abspath(f)
			files.append(dirc + '/' + f)
	return files
def GetOptions():
	parser = OptionParser()
	parser.add_option('-','--',dest = '', metavar = '', help = '')
	(options,args) = parser.parse_args()
	
	return
def GetBaseName(fname):
	return '.'.join(fname.split('.')[:-1])
def main():
	return	

if __name__=='__main__':
	main()
