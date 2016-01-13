"""



"""

import sys
import os

import auxiliary as aux

import inspect as ins

import ROOT as rt
import numpy as np
import datetime
import calendar


#gROOT.Reset()
gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
gStyle.SetOptFit(111111)
gStyle.SetPalette(1)



class PlotObject(object):
	pass


class Dataset(object):
	pass


class Histo(object):
	'''
	for histogram objects
	'''
	def __init__(self,TH):
		self.name = strName
		
		
	def __str__(self):
		return self.name 



class HistoCollection(object):
	'''
	collects histograms
	'''
	def __init__(self):
		self.count = 0
		self.h = {}
	
	
	def add(self,TH):
		self.count += 1
		
		
	
	def plotlist(self):
		# estimate the plot order of the histograms
		# using their integral and returns a list with the correct order
		pass
	





class myCanvas(object):
	'''
	create a canvas
	
	'''
	def __init__(self,name):
		self.opt = 0
		
		self.rat = 0.3
		
		self.cwidth = 800
		self.cheight = 900
		
		self.cname = name
		
		# ranges for the x axis
		self.xmin = -1
		self.xmax = 0
		
		#self.legend = TLegend(0.9,0.7,0.65,0.9)
		
	
	#def create(self):
	#	return TCanvas("name","name",800,900)
	
	def createCanvas(self, legend):
		self.l = legend
		
		self.htem = self.data.Clone("htem")
		self.htem.SetMinimum(-(self.htem.GetMaximum()/10.))
		
		self.createRatio()
		
		#self.f1 = self.ratio.getOneline()
		if(self.xmin != -1 and self.xmax != 0):
			self.mc.GetXaxis().SetRange(self.xmin, self.xmax)
			self.data.GetXaxis().SetRange(self.xmin, self.xmax)
			self.htem.GetXaxis().SetRange(self.xmin, self.xmax)
			#self.stack.GetXaxis().SetRange(self.xmin, self.xmax)
		
		self.c = TCanvas(self.cname,self.cname,self.cwidth,self.cheight)
		self.pad1 = TPad("pad1","pad1",0,self.rat,1,1)
		self.pad2 = TPad("pad2","pad2",0,0,1,self.rat)
		self.pad1.SetBottomMargin(0)
		self.pad2.SetTopMargin(0)
		
		self.pad1.Draw()
		self.pad1.cd()
		
		
		
		if(self.data.GetMaximum() > self.mc.GetMaximum()):
			print "data>mc"
			#self.htem = self.data.Clone()
			#self.htem.SetMinimum(-(self.htem.GetMaximum()/10.))
			self.htem.Draw("ep")
			
			
			#self.data.SetMinimum(-(self.data.GetMaximum()/10.))
			#self.data.Draw("ep")
			self.stack.Draw("same hist")
			self.data.Draw("same ep")
		else:
			print "mc>data"
			self.mc.SetMinimum(-(self.mc.GetMaximum()/10.))
			self.mc.Draw("hist")
			self.stack.Draw("same hist")
			self.data.Draw("same ep")
			
		self.l.Draw("same")
		
		self.c.cd()
		
		self.pad2.Draw()
		self.pad2.cd()
		self.ratio.getHis().Draw("ep")
		self.ratio.getOneline().Draw("same")
		
		self.c.cd()
	
	def addMC(self,h):
		# needs a THStack
		# "last" is the stacked histogram
		self.stack = h
		self.mc = h.GetStack().Last()
		
	
	
	def addData(self,h):
		# needs a histogram
		self.data = h
		
		self.data.SetMarkerStyle(20)
	
	
	def createRatio(self):
		# ratio of MC/Data
		self.ratio = Ratio(self.mc,self.data)
		
		if(self.xmin != -1 and self.xmax != 0):
			self.ratio.xmin = self.xmin
			self.ratio.xmax = self.xmax
		
		
	def plotOrder(self):
		n1 = self.data.GetMaximum()
		n2 = self.mc.GetMaximum()
		if(n1>n2):
			return n1, n2 # first draw data
		else:
			return n2, n1 # first draw mc
		
	
	def Styles(self):
		pass
	
	
	def getCanvas(self):
		return self.c
	
	
	def savePlot(self):
		
		
		
		
		return 0
	


class myCanvasRatio(myCanvas):
	'''
	OLD
	special canvas with ration plot
	'''
	def __init__(self,THhis,THratio, opt):
		self.his = THhis
		self.ratio = THratio
		myCanvas.__init__(self,opt)
		
	
	
	def pads(self, rat):
		
		self.pad1 = TPad("pad1","pad1",0,rat,1,1)
		self.pad2 = TPad("pad2","pad2",0,0,1,rat)
		self.pad1.SetBottomMargin(0)
		self.pad2.SetBottomMargin(0)
		
		
	
	
	def create(self):
		
		self._c = TCanvas("_c","_c",800,900)
		self.pads(0.3)
		
		self.pad1.Draw()
		self.pad1.cd()
		
		




class Ratio(object):
	'''
	ratio class.
	both histograms need to have the same number of bins
	and x axis range.
	'''
	def __init__(self,THnum,THden):
		self.num = THnum
		self.den = THden
		
		self.numtitle = self.num.GetXaxis().GetTitle()
		self.dentitle = self.num.GetXaxis().GetTitle()
		
		self.ytitle = self.numtitle + " / " + self.dentitle
		self.xtitle = ""
		
		self.c1 = 1
		self.c2 = 1
		
		self.ymin = 0.0 # standard range
		self.ymax = None
		
		self.__xmin = self.num.GetXaxis().GetXmin()
		self.__xmax = self.num.GetXaxis().GetXmax()
		
		self.xmin = self.__xmin
		self.xmax = self.__xmax
		
		self.nbins = self.num.GetNbinsX()
		self.name = self.num.GetName()
		#self.title = self.name + ";" + self.xtitle + ";" + self.ytitle
		self.title = ";" + self.numtitle + ";ratio"
		
		
		self.hisratio = 0
	
	
	def initRatio(self):
		
		self.hisratio = TH1F(self.name,self.title,self.nbins,self.__xmin,self.__xmax)
		
	
	def setOpt(self):
		
		# y range:
		self.ymax = self.hisratio.GetMaximum() + 0.2
		if(self.ymax < 1.5):
			self.ymax = 1.5
		self.hisratio.SetMinimum(self.ymin)
		self.hisratio.SetMaximum(self.ymax)
		
		
		#self.hisratio.GetXaxis().SetRange(self.xmin, self.xmax)

		#self.hisratio.GetXaxis().SetLimits(self.xmin,self.xmax)
		# x range:
		
		self.hisratio.SetMarkerStyle(20) # full circle
	
	
	def division(self):
		# TH1F->Divide(TH1 *h1, TH1 *h2,c1 = 1,c2 = 1,Option_t *option = "")
		# this = c1*h1/(c2*h2)
		self.hisratio.Divide(	self.num,
								self.den,
								self.c1,
								self.c2)
	
	def rebin(self,c_):
		self.hisratio.Rebin(c_)
	
	def getHis(self):
		self.initRatio()
		self.division()
		self.setOpt()
		return self.hisratio
		
	def getOneline(self):
		self.f1 = TF1("f1","1.",self.__xmin,self.__xmax)
		self.f1.SetLineWidth(1)
		return self.f1
	
	def hisdraw(self):
		pass
	




class Eff(object):
	'''
	for efficiency objects
	'''
	def __init__(self, name):
		self.var = 0
		
		self.ymin = 0
		self.ymax = 0.05 # 
		
		self.cwidth = 800
		self.cheight = 800
		
		self.cname = name
		
		self.graph = 0
		self.ffit = 0
	
	def addPassed(self,h):
		self.hpassed = h
	
	def addTotal(self, h):
		self.htotal = h
	
	def Graph(self):
		if(TEfficiency.CheckConsistency(self.hpassed,
										self.htotal)):
			print("consistency checked")
			self.eff = TEfficiency(	self.hpassed,
									self.htotal)
			
			self.graph = self.eff.CreateGraph()
			
			#return self.graph
		else:
			print("not consistent")
	
	def getGraph(self):
		if self.graph == 0:
			self.Graph()
		
		return self.graph
	
	def fit(self,func):
		self.ffit = func
		if self.graph == 0:
			self.Graph()
		
		self.graph.Fit(self.ffit,"r")
		
	
	
	def Canvas(self):
		'''
		 get a canvas with the efficiency
		'''
		if self.graph == 0:
			self.Graph()
		
		self.graph.SetTitle(self.cname)
		self.graph.SetMinimum(self.ymin)
		self.graph.SetMaximum(self.ymax)
		
		self.c = TCanvas(self.cname,self.cname,self.cwidth,self.cheight)
		self.c.cd()
		self.graph.Draw("ap")
		if self.ffit != 0:
			self.ffit.Draw("same")
		
		
		
	
	






class GraphText(object):
	def __init__(self,lum):
		self.lum = lum
	





def Fakerate(object):
	def __init__(self,THpassed,THtotal):
		self.passed = THpassed
		self.total = THtotal
		self.eff = Eff(THpassed,THtotal)
	
	def plot(self, c):
		return 0
	
def GiveTimestamp(strIn):
	'''
	returns the timestamp of a given string,
	if the structure of the string is of the form:
		"... _xyz_timestamp.root"
	'''
	return True
	
def GiveOutputString(strIn):
	'''
	takes the timestamp from the given string to
	create a string for the matching plot output file.
	an additional timestamp 144XXXXXX is added.
	
	input string hast to be of the form:
		../selector_rootFiles/my_selector_results_DY_v1_1447935068.root
		
	return:
		../plot_rootFiles/2015_48/my_plot_results_v1_1447935068_144XXXXXX.root
	'''
	ss = "../plot_rootFiles/"
	foldername = GiveYearWeekFoldername()
	strTemp = strIn.split("_")	
	strPath = ss + foldername + "my_plot_results_"
	dd = datetime.datetime.utcnow()
	strTimestamp = str(calendar.timegm(dd.timetuple()))
	strFormerTitle = strTemp[len(strTemp)-2]
	strEnd = strTemp[len(strTemp)-1] # 14XXXXXXXXXXX.root
	strFormerTimestamp = strEnd.split(".")[0]
	strFile = strFormerTitle + "_" + strFormerTimestamp + "_" + strTimestamp + ".root"
	if ensure_Dir(ss+foldername):
		#print "dir ensured"
		return strPath + strFile
	else:
		#print "dir not ensured"
		return False


def inputFileNumberOfLines():
	f = "outputfiles.txt"
	k = 0
	with open(f,"r") as s:
		for line in s:
			k += 1
	return k
	

def inputFile(n=""):
	'''
	returns the last line (default) or given line of the file outputfiles.txt,
	cleaned by the \n newline sign
	'''
	f = "outputfiles.txt"
	k = 0
	l = False
	if(n!=""):
		l = True
	with open(f,"r") as s:
		for line in s:
			k += 1
			if(l and k == n):
				return line.replace("\n","")
	return line.replace("\n","")


def FilelistDatetime():
	'''
	reads the file outputfiles.txt used by plot_tool.py
	and writes a list with the filenames and the 
	datetime given from the timestamp, also the file size in byte
	
	output: files_date.txt
	'''
	fi = "outputfiles.txt"
	fo = "files_date.txt"
	lines = []
	with open(fi,"r") as s:
		for line in s:
			path = line.replace("\n","")
			try:
				fsize = str(os.path.getsize(path))
			except os.error:
				fsize = " n/a "
			ss = line.split("_")
			d = int(ss[len(ss)-1].split(".")[0])		#date from filename
			ds = str(datetime.datetime.fromtimestamp(d)) #datestring
			lines.append(ds + "\t" + fsize + "\t" + line)
	with open(fo,"w") as s:
		for l in lines:
			s.write(l)
	return True


def GivePlotfolder():
	'''
	returns the path of the current plot folder (year and week number)
	example: "../plots/2015_48/" (including the last / )
	'''
	ss = "../plots/"
	foldername = GiveYearWeekFoldername()
	if ensure_Dir(ss+foldername):
		#print "dir ensured"
		return ss+foldername
	else:
		#print "dir not ensured"
		return False

def GiveYearWeekFoldername():
	'''
	returns a string constructed by year and week, e.g.:
	"2015_48/"
	'''
	today = datetime.date.today()
	y,w,d = today.isocalendar()
	return str(y)+"_"+str(w)+"/"

def GivePlotFilePath(strFile):
	'''
	returns the file
	'''
	path = GivePlotfolder()
	return path+strFile


def Folder():
	'''
	
	'''
	pass

def ensure_Dir(path):
	'''
	checks if a path exists, otherwise create it
	'''
	d = os.path.dirname(path)
	if not os.path.exists(d):
		# print "path doesnt exist"
		os.mkdir(d)
		return True
	elif os.path.exists(d):
		# print "path exists"
		return True


def someTest():
	return 0






def give_hCutFlow(path):
	'''
	gives back the first bin of hCutFlow
	'''
	data = TFile(path)
	hcutflow = data.Get("TreeWriter/hCutFlow")
	return hcutflow.GetBinContent(2)



def giveDataPath(key):
	'''
	
	'''
	datasets = [	"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",
					"TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_v01.root",
					"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",
					"ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",
					"DoubleEG_Run2015D-05Oct2015-v1.root",
					"DoubleEG_Run2015D-PromptReco-v4.root"
				]
	path = "/user/rmeyer/mergedNTuples/";#
	for k in datasets:
		if key in k:
			return path+k
	return False


























