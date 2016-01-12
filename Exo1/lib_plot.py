"""
some classes and functions

"""

import ROOT as rt
from lib_filemanagement import *
import numpy as np
import datetime
import calendar
import sys
import os
from array import array


#gROOT.Reset()
rt.gROOT.SetStyle("Plain")
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(111111)
rt.gStyle.SetPalette(1)



class FileSaver():
	'''
	
	'''
	def __init__(self, filename, treepath, *args):
		# args:
		# 0: text plotted
		# 1:
		 
		self.filename = filename # 
		self.treepath = treepath # e.g. 'HtBins'
		
		
	
	
	def save(self, obj, name):
		name = self.treepath+'/'+name
		if isinstance(obj, rt.TCanvas):
			obj.Update()
			
		
		f = rt.TFile(self.filename, "update")
		f.cd()
		
		
		
		
		
	




class PlotObject(object):
	pass


class Dataset(object):
	pass


class Histo(object):
	'''
	for histogram objects
	'''
	def __init__(self):
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
	
	
	def add(self,h):
		self.count += 1
		
		
	
	def plotlist(self):
		# estimate the plot order of the histograms
		# using their integral and returns a list with the correct order
		pass
	



class ratioCanvasHF(object):
	'''
	create a canvas with a ratio plot
	from a histogram to a fitted function
	where h=data, f=fit, ratioplot=(h-f)/sigma_h
	
	standard size: width*height = 
	
	usage:
	cr = ratioCanvasFunction(h,f, "name")
	cr.draw() # create the canvas
	
	'''
	def __init__(self,h,f,name,csize=900):
		
		self.cname = name
		
		self.cwidth = csize
		self.cheight = csize-csize/4
		
		self.h = h
		self.f = f
		
		self.h.SetMarkerColor(rt.kRed)
		self.h.SetLineWidth(1)
		self.h.SetMarkerStyle(20)
		
		self.xmin = h.GetXaxis().GetXmin()
		self.xmax = h.GetXaxis().GetXmax()
		
		self.h.SetTitle(name)
		
		# for graphic reasons, overwritten when createRatio is called
		self.bin0 = 0
		self.binEnd = 0
		
		# create ratio plot histogram
		self.ratio = self.createRatio(h, f)
		
		# percentage of the ratio plot
		self.rat = 0.37
		
		
	
	def draw(self, *log):
		
		if log:
			pass
		else:
			log = ["logx","logy"]
		
		
		self.c = rt.TCanvas(self.cname,self.cname,self.cwidth,self.cheight)
		self.c.cd()
		
		if "EBEB" in self.cname:
			#print "EBEB"
			img = rt.TImage.Open("CMS-PAS-EXO-15-004_Figure_003-a.png")
			#img.Draw("x")
		
		if "EBEE" in self.cname:
			#print "EBEE"
			img = rt.TImage.Open("CMS-PAS-EXO-15-004_Figure_003-b.png")
			#img.Draw("x")
		
		
		# zero line
		self.f0 = rt.TF1("f1","0",self.xmin,self.xmax)
		self.f0.SetLineWidth(2)
		
		#self.c.SetBottomMargin(0.1)
		
		# 
		self.pad1 = rt.TPad("pad1","pad1",0,self.rat,1,1)
		self.pad2 = rt.TPad("pad2","pad2",0,0.0,1,self.rat)
		
		self.pad1.SetBottomMargin(0.01)
		
		self.pad2.SetTopMargin(0)
		self.pad2.SetBottomMargin(0.15)
		
		self.pad1.SetFillStyle(4000)
		self.pad1.SetFrameFillStyle(4000)
		self.pad2.SetFillStyle(4000)
		self.pad2.SetFrameFillStyle(4000)
		
		self.pad1.SetGridx()
		self.pad1.SetGridy()
		
		self.pad2.SetGridx()
		self.pad2.SetGridy()
		
		
		
		if "logx" in log:
			self.pad1.SetLogx()
			self.pad2.SetLogx()
		
		if "logy" in log:
			self.pad1.SetLogy()
			self.h.SetMinimum(10e-3)
			#self.h.SetMinimum(self.f.Eval(1.01*self.binEnd))
			self.h.SetMaximum(300)
			#self.h.SetMaximum(self.f.Eval(0.5*self.bin0))
			
		self.c.cd()
		
		self.h.SetLineWidth(2)
		
		#self.h.SetLineStyle(8)
		
		# pad1: histogram + fit
		self.pad1.Draw()
		self.pad1.cd()
		self.h.DrawCopy("hist")
		self.h.SetLineWidth(1)
		self.h.Draw("ep same")
		self.f.Draw("same")
		
		self.c.cd()
		
		# pad2: ratioplot (+ zero line)
		self.pad2.Draw()
		self.pad2.cd()
		self.ratio.Draw("ep")
		self.f0.Draw("same")
		self.c.cd()
		
		
		
		
		#self.c.SaveAs("plots/" + self.cname + ".png")
		#self.c.SaveAs("plots/" + self.cname + ".root")
		
		# append plot canvases to root tree
		filename = "plots/plots.root"
		myfile = rt.TFile(filename,"UPDATE")
		self.c.Write()
		
		# '''
		
		#return self.c
	
	def addHisto(self,his, hf):
		csize = 900
		
		entries = his.GetNbinsX()
		
		xmin = his.GetXaxis().GetXmin()
		xmax = his.GetXaxis().GetXmax()
		
		binwidth = int((xmax-xmin)/entries)
		
		hname = his.GetName() +"_"+str(binwidth)
		
		#self.h.SetMaximum(h.GetMaximum())
		c = rt.TCanvas(hname,hname,csize,csize*3/4)
		c.SetLogx()
		c.SetLogy()
		#h.SetLineColor(kBlue)
		#self.h.Draw("ep")
		
		his.SetLineColor(rt.kBlack)
		his.SetLineWidth(2)
		hf.SetLineColor(rt.kBlue)
		hf.SetLineWidth(2)
		
		c.cd()
		his.Draw("hist")
		hf.Draw("same hist")
		self.h.Draw("same ep")
		c.Update()
		
		
		
		# append plot canvas to root tree
		filename = "plots/cumulatives.root"
		myfile = rt.TFile(filename,"UPDATE")
		c.Write()
		#c.SaveAs("plots/firstone.root")
		
		# '''
		
	
	def getHist(self):
		return self.h
	
	def createRatio(self,h,f):
		# creates a ratio histogram Diff = (Data-Fit)/Sigma_stat
		# from "data" histogram(h) and fit(f)
		#
		# WARNING: the errors of Diff are estimated gaussian:
		# err = dDiff/dData*Sigma_stat = 1/Sigma_stat*Sigma_stat = 1	!!!
		# 
		entries = h.GetNbinsX()
		
		diff =[]
		err = []
		xvals = []
		for i in range(0,entries):
			data = h.GetBinContent(i)
			terr = np.sqrt(data)
			if terr == 0.: # to avoid 1/NULL or errors
				terr = 1.
			err.append(terr)
			x = h.GetBinCenter(i)
			xvals.append(x)
			fitted = f.Eval(x)
			diff.append((data-fitted)/terr)
		
		self.bin0 = xvals[0]
		self.binEnd = xvals[entries-1]
		
		htitle = ";" + h.GetXaxis().GetTitle() + ";(data-fit)/#sigma_{stat}"
		hname = "ratio" + self.cname
		ratio = rt.TH1F(hname,htitle,entries,self.xmin,self.xmax)
		for i in range(0, entries):
			ratio.SetBinContent(i, diff[i])
			ratio.SetBinError(i,1.)
		
		ratio.GetXaxis().SetLabelSize(0.07)
		ratio.GetXaxis().SetTitleSize(0.07)
		ratio.GetYaxis().SetLabelSize(0.05)
		ratio.GetYaxis().SetTitleSize(0.04)
		
		ratio.SetMarkerStyle(20)
		ratio.SetMarkerColor(rt.kRed)
		
		ratio.SetStats(0) # dont plot the stat box from the ratio
		
		return ratio




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
		
		self.c = rt.TCanvas(self.cname,self.cname,self.cwidth,self.cheight)
		self.pad1 = rt.TPad("pad1","pad1",0,self.rat,1,1)
		self.pad2 = rt.TPad("pad2","pad2",0,0,1,self.rat)
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
		
		self.pad1 = rt.TPad("pad1","pad1",0,rat,1,1)
		self.pad2 = rt.TPad("pad2","pad2",0,0,1,rat)
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
		
		self.hisratio = rt.TH1F(self.name,self.title,self.nbins,self.__xmin,self.__xmax)
		
	
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
		self.f1 = rt.TF1("f1","1.",self.__xmin,self.__xmax)
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
		if(rt.TEfficiency.CheckConsistency(self.hpassed,
										self.htotal)):
			print("consistency checked")
			self.eff = rt.TEfficiency(	self.hpassed,
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
		
		self.c = rt.TCanvas(self.cname,self.cname,self.cwidth,self.cheight)
		self.c.cd()
		self.graph.Draw("ap")
		if self.ffit != 0:
			self.ffit.Draw("same")
		
		
		
	
	






class GraphText(object):
	def __init__(self,lum):
		self.lum = lum
	





class Fakerate(object):
	def __init__(self,THpassed,THtotal):
		self.passed = THpassed
		self.total = THtotal
		self.eff = Eff(THpassed,THtotal)
	
	def plot(self, c):
		return 0





def createHistoFromFunction(f, h):
	# create a histogram from a given function 
	# 
	fname = f.GetName()
	entries = h.GetNbinsX()
	hname = h.GetName()
	xmin = h.GetXaxis().GetXmin()
	xmax = h.GetXaxis().GetXmax()
	xtitle = h.GetXaxis().GetTitle()
	ytitle = h.GetYaxis().GetTitle() 
	binwidth = int((xmax-xmin)/entries)
	
	hf = rt.TH1F(hname+"_"+str(binwidth),"cumulative "+fname+";"+xtitle+";"+ytitle,entries,xmin,xmax) 
	
	# fill histogram with the function values at bin centers
	for i in range(0,entries):
		x = h.GetBinCenter(i)
		hf.SetBinContent(i, f.Eval(x))
	
	return hf


def CumulativePlot(h, direction):
	# creates a cumulative histogram
	# direction:	"up" from lower x add up bins to upper x
	#				"down" from upper x add up bins to lower x
	entries = h.GetNbinsX()
	hname = h.GetName()
	xmin = h.GetXaxis().GetXmin()
	xmax = h.GetXaxis().GetXmax()
	xtitle = h.GetXaxis().GetTitle()
	ytitle = h.GetYaxis().GetTitle()
	
	binwidth = int((xmax-xmin)/entries)
	
	summed = 0
	
	hc = rt.TH1F(hname+"_"+str(binwidth),"cumulative "+direction+";"+xtitle+";"+ytitle,entries,xmin,xmax)
	
	if direction == "up":
		r = range(0, entries-1, 1)
		summed = h.GetBinContent(1)
	if direction == "down":
		r = range(entries-1, 0, -1)
		summed = h.GetBinContent(entries)
	else:
		print "... please give correct direction"
		sys.exit()
	
	for i in r:
		summed = summed + h.GetBinContent(i)
		hc.SetBinContent(i, summed)
	
	return hc





def give_hCutFlow(path,bin_=2):
	'''
	gives back the first bin of hCutFlow
	'''
	data = rt.TFile(path)
	hcutflow = data.Get("TreeWriter/hCutFlow")
	return hcutflow.GetBinContent(bin_)



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


def giveMClumiScalefactor(key):
	'''
	returns the lumi scale factor for a monte carlo dataset
	'''
	datapath = giveDataPath(key)
	N = give_hCutFlow(datapath)
	sigma = giveCrossSectionMC(key)
	L = giveLumiData(key)
	return float(L)/float(N)*float(sigma)
	


def giveCrossSectionMC(key):
	'''
	returns cross section in pb
	'''
	cs = {		"DYJetsToLL": 6025.2, # NLO
				# "DYJetsToLL": 6104., # LO
				"TTGJets": 3.697, # NLO
				"WGToLNuG": 489., # NLO
				"ZGTo2LG": 117.864 # NLO
			}
	for k in cs.keys():
		if key in k:
			return cs[k]
		
	return False


def giveLumiData(key):
	'''
	returns lumi in pb-1
	'''
	lumi = {	"Run2015-05Oct2015": 2.11e3,		# originally in fb-1
				"Run2015-PromptReco": 2.11e3		# originally in fb-1
			}
	#or k in lumi.keys():
	#	if key in k:
	#		return lumi[k]
	val = 2.11e3
	return val


































































