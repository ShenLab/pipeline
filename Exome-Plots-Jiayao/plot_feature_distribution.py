#!/home/local/users/jw/bin/python2.7
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# Plot features in a features table file for consensus analysis.
#========================================================================================================

from optparse import OptionParser
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.ticker as ticker
import math
from scipy.stats import gaussian_kde

def AdHoc_SOR(A,B,C,D):
	A,B,C,D = map(float,[A,B,C,D])
	try:
		R = (A*D)/(B*C)
		return R + 1/R
	except ZeroDivisionError:
		return 0

def GetOptions():
	parser = OptionParser()
	parser.add_option('-i','--input',dest = 'INPUT', metavar = 'INPUT', help = 'input features table')
	parser.add_option('-f','--feature',dest = 'FEATURE', metavar = 'FEATURE', help = 'input features to be plot')
	(options,args) = parser.parse_args()
	return options.INPUT,options.FEATURE

#========================================================================================================
# Plot Frequency Distribution
#========================================================================================================
def pick_PL(PL):
	PL = map(float,PL.split(','))
	PL.sort()
	return PL[1]
def pick_GL(GL):
	GL = map(float,GL.split(','))
	GL.sort()
	return GL[-2]

def Get_data_freq(InpFil,feature):
	fin = open(InpFil,'rb')
	headers = fin.readline().strip().split('\t')
	if feature == 'gatk-fmt-DP':
		index = headers.index(feature)
		data = []
		for l in fin:
			llist = l.strip().split('\t')
			tmp = llist[index]
			if tmp != '.':
				tmp = float(tmp)
			#if tmp < 20:
				data.append(tmp)
	elif feature == 'st-fmt-AD':
		index = headers.index(feature)
		data = []
		for l in fin:
			llist = l.strip().split('\t')
			try:
				tmp1,tmp2 = map(float,llist[index].split(','))
				tmp = tmp2/(tmp1+tmp2)
				if tmp != 1:
					data.append(tmp)
			except:
				#print llist[index].split('\t')
				continue
	elif feature == 'st-fmt-PL':
		index = headers.index(feature)
		data = []
		for l in fin:
			llist = l.strip().split('\t')
			try:
				tmp = pick_PL(llist[index])
				data.append(tmp)
			except:
				#print llist[index].split(',')
				continue
	else:
		index = headers.index(feature)
		data = []
		for l in fin:
			llist = l.strip().split('\t')
			try:
				tmp = float(llist[index])
				data.append(tmp)
			except:
				#print llist[index].split('\t')
				continue
	return data

def Plot_one_freq(InpFils,feature):
	data_freq_1 = Get_data_freq(InpFils[0],feature)
	data_freq_2 = Get_data_freq(InpFils[1],feature)
	data_freq_3 = Get_data_freq(InpFils[2],feature)
	data_freq_4 = Get_data_freq(InpFils[3],feature)
	data_freq = [data_freq_1,data_freq_2,data_freq_3,data_freq_4]
	plot_freq(data_freq,feature)
def plot_freq(datas,feature):
	colors = ['red','yellow','blue','black']
	for i,(data) in enumerate(datas):
		#print data
		density = gaussian_kde(data)
		xs = np.linspace(-40,40,100)
		density.covariance_factor = lambda : .10
		density._compute_covariance()
		plt.plot(xs,density(xs),linewidth=2.0)
		plt.legend(['1/4', '2/4', '3/4', '4/4'], loc='upper right')
	#plt.xticks(np.arange(0, 1, 0.1),rotation='vertical')
	plt.xlabel('gatk-VQSLOD')
	plt.ylabel('Frequency')
	plt.title(feature+' For All Variants')
	#plt.title(feature+' For Rare Variants')
	plt.show()
#========================================================================================================
# Plot Scatter Plot (For NR-NV;)
#========================================================================================================
def Plot_one_scatter(InpFils,feature):
	x_1,y_1 = Get_data_scatter(InpFils[0],feature)
	x_2,y_2 = Get_data_scatter(InpFils[1],feature)
	x_3,y_3 = Get_data_scatter(InpFils[2],feature)
	x_4,y_4 = Get_data_scatter(InpFils[3],feature)
	xx = [x_1,x_2,x_3,x_4]
	yy = [y_1,y_2,y_3,y_4]
	#plot_scatter(xx,yy,feature)
	plot_scatter_subs(xx,yy,feature)
def Get_data_scatter(InpFil,feature):
	fin = open(InpFil,'rb')
	headers = fin.readline().strip().split('\t')
	index = headers.index(feature)

	index_1 = headers.index('fb-SRP') # For Platypus
	index_2 = headers.index('fb-SAP') # For Platypus

	#index_1 = headers.index(feature)
	#index_2 = headers.index('ExAC_ALL')

	x,y = [],[]
	for l in fin:
		llist = l.strip().split('\t')
		"""
		tmp = llist[index]
		if tmp != '.':
			tmp = tmp.split(',')
			#if 1:
			if int(tmp[0]) <= 150 and int(tmp[1]) <= 150:
				x.append(float(tmp[0])) #GATK
				#x.append(float(tmp[0])-float(tmp[1]))
				y.append(float(tmp[1]))
		"""
		tmp_1 = llist[index_1].split(',')[0]
		tmp_2 = llist[index_2].split(',')[0]
		#if tmp_1 != '.' and tmp_2 != '.' and int(tmp_1) <= 200:
		if tmp_1 != '.':

			tmp_2 = 0 if tmp_2 == '.' else tmp_2
			
			if float(tmp_2) <= 150:
				#x.append(float(tmp_1)-float(tmp_2))
				x.append(float(tmp_1))
				y.append(float(tmp_2))
		
	return x,y
def plot_scatter(xx,yy,feature):
	makers = ['+','+','+','+']
	colors = ['red','yellow','blue','green']
	labels = ['1/4', '2/4', '3/4', '4/4']
	for i,(x,y) in enumerate(zip(xx,yy)):

		plt.scatter(x,y, color = colors[i], s=1, label=labels[i], marker='o') 
	plt.title(feature)
	plt.xlabel('fb-SRP')
	plt.ylabel('fb-SAP')

	#plt.xscale('log')
	#plt.yscale('log')
	#plt.xticks(np.arange(min(x), max(x)+1),0.0001) 
	#plt.yticks(np.arange(min(x), max(x)+1),0.0001)
	plt.grid(True)
	plt.legend(loc = 'upper right') 
	plt.show()
def plot_scatter_subs(xx,yy,feature):
	makers = ['+','+','+','+']
	colors = ['red','grey','blue','green']
	labels = ['1/4', '2/4', '3/4', '4/4']
	subs = [221,222,223,224]
	for i,(x,y) in enumerate(zip(xx,yy)):
		#x = [math.log10(i) for i in x]
		#y = [math.log10(i) for i in y]
		plt.subplot(subs[i])
		plt.scatter(x,y, color = colors[i], s=2, label=labels[i], marker='x') 
		#plt.xlabel(feature)
		#plt.ylabel('ExAC_ALL')

		plt.xlabel('fb-SRP')
		plt.ylabel('fb-SAP')
		plt.title(feature+': '+labels[i])


		axes = plt.gca()
		#axes.set_xlim([0,150])
		#axes.set_ylim([0,150])
		plt.margins(0.2)
		#plt.xscale('log')
		#plt.yscale('log')
		#plt.xticks(np.arange(0, max(x), 0.05),rotation='vertical')
		#plt.xticks(np.arange(0, max(x), 0.1))
		#plt.yticks(np.arange(0, max(y), 0.001))
		plt.grid(True)

	#plt.legend(loc = 'upper right') 
	plt.show()

#========================================================================================================
# Plot Scatter Plot (For )
#========================================================================================================
def Plot(InpFils,Features):
	for feature in Features:
		Plot_one_freq(InpFils,feature)
		#Plot_one_scatter(InpFils,feature)

def main():
	#InpFil,Feature = GetOptions()
	InpFils,Feature = GetOptions_2()
	Plot(InpFils,Feature)

def GetOptions_2():

	one = 'uniqall_features.tsv'
	two = '2outof4all_features.tsv'
	three = '3outof4all_features.tsv'
	four = '4outof4all_features.tsv'


	indel = './indel/'
	snv = './snv/'
	#indel = './rare/indel/'
	#snv = './rare/snv/'
	files = [snv + tmp for tmp in [one,two,three,four]]
	#files = [indel + tmp for tmp in [one,two,three,four]]

	#features = ['gatk-fmt-DP','st-fmt-DP','pt-fmt-NR','fb-fmt-DP']
	#features = ['fb-fmt-GL']
	#features = ['gatk-MLEAF','st-AF1','pt-FR','fb-AF']
	features = ['gatk-VQSLOD']
	return files,features

if __name__=='__main__':
	main()

