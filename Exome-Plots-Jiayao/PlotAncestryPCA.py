# coding: utf-8
import argparse
import csv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random
from mpl_toolkits.mplot3d import Axes3D

ControlDefault="/home/local/users/jw/resources/AncestryPCA/resources/1KG_AJ_Domi_PCAcontrol.vcf.gz"
IndvPanelDefault="/home/local/users/jw/resources/AncestryPCA/resources/AncestryPCA.master.panel"

class Control:
	def __init__(self, row):
		self.ID, self.Pop, self.SuperPop, self.Gender = row

class Evec:
	def __init__(self, row):
		self.ID = row[0].split(":")[-1]
		self.Partition = row[-1]
		self.EigenValues = row[1:-1]

class PloyAncestryPCA:
	def __init__(self, ControlPanel, SampleEVEC, SampleEVAL):
		self.Name = SampleEVEC.split("/")[-1].split(".")[0]
		self.ControlPanelFil = ControlPanel
		self.EvecFil = SampleEVEC
		self.EvalFil = SampleEVAL

	def run(self):
		self.loadControlPanel()
		self.loadEigenEval()
		self.loadEigenValue()
		self.Plot(self.Name)
		self.Plot3D(self.Name)

	def loadControlPanel(self):
		# sample  pop     super_pop       gender
		self.ControlPanel = {}
		fin = open(self.ControlPanelFil, 'rb')
		reader = csv.reader(fin, delimiter='\t')
		header = reader.next()
		for row in reader:
			Indv = Control(row)
			self.ControlPanel[Indv.ID] = Indv

	def loadEigenEval(self):
		fin = open(self.EvalFil)
		self.EigenEvals = []
		SUM = 0
		for l in fin:
			tmp = float(l.strip().split()[0])
			SUM += tmp
			self.EigenEvals.append(tmp)
		for i,tmp in enumerate(self.EigenEvals):
			self.EigenEvals[i] = tmp/SUM

	def loadEigenValue(self):
		fin = open(self.EvecFil, 'rb')
		self.Colors = {'SAS':"green" , 'EAS':"blue", 'AMR':"orange", 'AFR':"purple", 'EUR':"red", 'AJ':"
				" 'CASE':"grey"}
		self.SuperPops = {'SAS':[] , 'EAS':[], 'AMR':[], 'AFR':[], 'EUR':[], 'AJ':[], 'DOMI':[]} #'SAS', 'EAS', 'AMR', 'AFR', 'EUR'
		self.Cases = []
		header = fin.next()
		for line in fin:
			row = line.strip().split()
			evec = Evec(row)
			if evec.Partition == "Case":
				self.Cases.append(evec)
			elif evec.Partition == "Control":
				try:
					superpop = self.ControlPanel[evec.ID].SuperPop
					self.SuperPops[superpop].append(evec)
				except KeyError:
					print "Sample %s marked as control but not in Control Panel" % (evec.ID)


	def ReduceSuperPop(self):
		self.Colors = {'SAS':"green" , 'EAS':"blue", 'AMR':"orange", 'AFR':"purple", 'EUR':"red", 'CASE':"grey"}
		self.SuperPops = {'SAS':[] , 'EAS':[], 'AMR':[], 'AFR':[], 'EUR':[]} #'SAS', 'EAS', 'AMR', 'AFR', 'EUR'
		for k,v in self.ControlPanel.items():
			superpop = v.SuperPop
			if superpop not in self.SuperPops:
				print "{}'s SuperPop {} is not included in default SuperPop: {}".format(k, superpop, "SAS, EAS, AMR, AFR, EUR")
			else:
				self.SuperPops[superpop].append(k)

	def ReadEigenvaluesFromList(self, List, idx1, idx2):
		X, Y = [], []
		for evec in List:
			X.append(evec.EigenValues[idx1])
			Y.append(evec.EigenValues[idx2])
		return X, Y

	def Plot(self, OutName):
		with PdfPages('%s.pdf'%OutName, 'wb') as pdf:
			for PCx, PCy in [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10)]:
				plt.figure(figsize=(15, 15))
				for k,SuperPop in self.SuperPops.items():
					color = self.Colors[k]
					X, Y = self.ReadEigenvaluesFromList(SuperPop, PCx-1, PCy-1)
					plt.scatter(X, Y, color=color, marker='x', label=k)
				X, Y = self.ReadEigenvaluesFromList(self.Cases, PCx-1, PCy-1)
				plt.scatter(X, Y, color="grey", marker='o', label="Cases")
				plt.xlabel("PC%d, %.3f%s of variance" % (PCx,self.EigenEvals[PCx-1]*100,"%"))
				plt.ylabel("PC%d, %.3f%s of variance" % (PCy,self.EigenEvals[PCy-1]*100,"%"))
				plt.title('PC%d - PC%d'%(PCx,PCy))
				plt.legend(loc='upper right')
				plt.grid(True)
				#plt.show()
				pdf.savefig()
				plt.close()

	def Plot3D(self, OutName):
		with PdfPages('%s.3D.pdf'%OutName, 'wb') as pdf:
			fig = plt.figure(figsize=(20, 20))
			ax = fig.add_subplot(111, projection='3d')
			# Plot Cases
			X, Y, Z = [], [], []
			for evec in self.Cases:
				X.append(float(evec.EigenValues[0]))
				Y.append(float(evec.EigenValues[1]))
				Z.append(float(evec.EigenValues[2]))
			ax.scatter(X, Y, Z, c="grey", marker='o', label="Cases")
			# Plot Controls
			for k,SuperPop in self.SuperPops.items():
				color = self.Colors[k]
				X, Y, Z = [], [], []
				for evec in SuperPop:
					X.append(float(evec.EigenValues[0]))
					Y.append(float(evec.EigenValues[1]))
					Z.append(float(evec.EigenValues[2]))
				ax.scatter(X, Y, Z, c=color, marker='x',label=k)

			ax.set_xlabel('PC%d, %.3f%s of variance'%(1,self.EigenEvals[0]*100,"%"))
			ax.set_ylabel('PC%d, %.3f%s of variance'%(2,self.EigenEvals[1]*100,"%"))
			ax.set_zlabel('PC%d, %.3f%s of variance'%(3,self.EigenEvals[2]*100,"%"))
			plt.title(OutName)
			plt.legend(loc='upper right')
			#plt.grid(True)
			#plt.show()
			pdf.savefig()
			plt.close()


def GetOptions():
	#ControlPanel = "integrated_call_samples_v3.20130502.ALL.panel"
	#EVEC = "data/RGN.plus.HapMap.pca.evec"
	#EVAL = "data/RGN.plus.HapMap.eval"
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--prefix', type=str, help='Prefix of PCA results')
	parser.add_argument('-c', '--control', type=str, default=IndvPanelDefault, help='Control Panel')
	parser.add_argument('--evec', type=str, help='EigenValue Evaluation File')
	parser.add_argument('--eval', type=str, help='Sample EigenVector Files')
	args = parser.parse_args()
	if args.prefix == None and (args.evec == None or args.eval == None):
		print "Missing Input Data. Please see help"
		exit()
	if args.prefix != None:
		args.evec = "%s.%s"%(args.prefix, "pca.evec")
		args.eval = "%s.%s"%(args.prefix, "eval")
	return args.control, args.evec, args.eval

def main():
	ControlPanel, SampleEVEC, SampleEVAL = GetOptions()
	instance = PloyAncestryPCA(ControlPanel, SampleEVEC, SampleEVAL)
	instance.run()

if __name__ == '__main__':
	main()
