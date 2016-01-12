
#from ROOT import *
import ROOT as rt
from lib_plot import *
from lib_filemanagement import *
from lib_auxiliary import *
from math import *
from array import array



# global objects
csize = 800


def fitMassSpectrum(x, pars):
	# diphoton spectrum function
	return pow(x[0],pars[0]+pars[1]*log(x[0]))
	


def bw(x, pars):
	# breit wigner
	return pars[0]*1./(pow((pars[1]-x[0]),2) + pow((x[0]*pars[2]),2))
	

def fit_tool(h, f, c):
	#h.Scale(1./float(h.GetEntries()))
	h.Fit(f)
	c.cd()
	f.Draw("same")
	
	return True


def fake_dep(fr, var):
	pass


def getDiphotonData(lines=2, lineoff=0):
	# read out data
	datakeys = {	"DYJetsToLL": 0,
					"TTGJets": 2,
					"WGToLNuG": 4,
					"ZGTo2LG": 6,
					"Run2015D-05Oct2015": 0,
					"Run2015D-PromptReco": 1
				}
	defcolor = {	0: rt.kGreen}
	
	# get number of lines from file:
	total_lines = inputFileNumberOfLines()
	h = {} #histograms
	
	h2 = {} # 2d histograms
	
	# read the files:
	for ll in range(total_lines-lineoff-lines+1,total_lines-lineoff+1):
		# filename
		strInputFile = inputFile(ll)
		for key in datakeys:
			if(key in strInputFile):
				hname = key
		#print hname
		h[hname] = {}
		h2[hname] = {}
		try:
			data = rt.TFile(strInputFile)
		except:
			return False
		
		# loop all histograms and search for EBEB and EBEE
		for key in data.GetListOfKeys():
			kgn = key.GetName()
			# search for the datatag
			if "_EB" in kgn:
				#print type(data.Get(kgn))
				h[hname][kgn] = data.Get(kgn)
				h[hname][kgn].SetDirectory(0)
			if "-EB" in kgn:
				h2[hname][kgn] = data.Get(kgn)
				h2[hname][kgn].SetDirectory(0)
			
		
		data.Close()
	
	return h, h2



# //////////////////////////////////////////////////////////////////////////////////////////////////
# plot mass spectrum
def massSpectrum(lines=2, lineoff=0):
	print sys.argv
	nrebin = 20
	
	# read out data
	h, h2 = getDiphotonData(lines, lineoff)
	
	
	# 1d:
	hEBEB = rt.TH1F("EBEB",";m [GeV];counts",2100,0,2100)
	hEBEE = rt.TH1F("EBEE",";m [GeV];counts",2100,0,2100)
	
	# add up both
	for key in h.keys():
		for qkey in h[key].keys():
			if "_EBEB" in qkey:
				hEBEB.Add(h[key][qkey])
			if "_EBEE" in qkey:
				hEBEE.Add(h[key][qkey])
	
	# 1d
	print "start binning plots:"
	for i in [5,10,12,14,15,16,18,20,22,24,25,26,28,30,32,34,35,36,38,40,45,50]:
		print "binwidth: ", i
		#plotterDiphoton(hEBEB, hEBEE, i)
	
	
	# handle 2d histograms
	hEBEBmet = rt.TH2F("metEBEB_bins",";m [GeV];met< [GeV]",
						2100,0,2100,
						2000,0,2000);
	hEBEBht = rt.TH2F("htEBEB_bins",";m [GeV];Ht< [GeV]",
						2100,0,2100,
						2000,0,2000)
	hEBEEmet = rt.TH2F("metEBEE_bins",";m [GeV];met< [GeV]",
						2100,0,2100,
						2000,0,2000);
	hEBEEht = rt.TH2F("htEBEE_bins",";m [GeV];Ht< [GeV]",
						2100,0,2100,
						2000,0,2000)
	
	for key in h2.keys():
		for qkey in h2[key].keys():
			#print h2[key][qkey]
			if "met-EBEB" in qkey:
				#print "met added!"
				hEBEBmet.Add(h2[key][qkey])
			if "ht-EBEB" in qkey:
				hEBEBht.Add(h2[key][qkey])
			if "met-EBEE" in qkey:
				hEBEEmet.Add(h2[key][qkey])
			if "ht-EBEE" in qkey:
				hEBEEht.Add(h2[key][qkey])
	
	# projection to TH1:
	# ProjectionX (const char *name="_px", Int_t firstybin=0, Int_t lastybin=-1, Option_t *option="") 
	hEBEBmetPx = {} # collections for different met and ht cuts
	hEBEBhtPx = {} # -"-
	hEBEEmetPx = {} # -"-
	hEBEEhtPx = {} # -"-
	
	
	
	# met:
	print "start met plots: "
	for met in range(0,201, 20):
		print "met cut: ", met
		hEBEBmetPx[met] = hEBEBmet.ProjectionX("MetEBEB_"+str(met), 0, met)
		hEBEEmetPx[met] = hEBEEmet.ProjectionX("MetEBEE_"+str(met), 0, met)
	
	plotterDiphoton(hEBEBmetPx[met], hEBEEmetPx[met])
	
	# ht:
	print "start ht plots: "
	for ht in range(0, 1001, 50):
		print "ht cut: ", ht
		hEBEBhtPx[ht] = hEBEBht.ProjectionX("HtEBEB_"+str(ht), 0, ht)
		hEBEEhtPx[ht] = hEBEEht.ProjectionX("HtEBEE_"+str(ht), 0, ht)
	
	plotterDiphoton(hEBEBhtPx[ht], hEBEEhtPx[ht])
	
	#print "hEBEBmetPx[100] : ",hEBEBmetPx[100].GetName()
	#print "hEBEEmetPx[100] : ",hEBEEmetPx[100].GetName()
	
	# 2d projections:
	#for met in range(0,201, 20):
		
	
	#for ht in range(0, 1001, 50):
		
	
	
	
	print "... massSpectrum() done"
	raw_input()
	
	return 0


def plotterDiphoton(hEBEB, hEBEE, *args):
	# needs two histograms
	
	#print "hEBEB: ", hEBEB.GetName()
	#print "hEBEE: ", hEBEE.GetName()
	
		# draw and fit:
	maxvalue = 1990 # for both
	
	binEB = 20 # standard binning: events per 20 GeV
	
	if len(sys.argv) > 1:
		binEB = int(sys.argv[1])
	if len(args) > 0:
		binEB = args[0]
	
	binEE = binEB
	
	# EBEB
	fmin = 230
	
	fmax = maxValFromBinwidth(fmin,maxvalue,binEB)
	
	hEBEB = rebin(hEBEB, range(fmin,fmax+1,binEB))
	hEBEB.GetXaxis().SetRangeUser(fmin, fmax)
	
	hEBEBc = CumulativePlot(hEBEB, "down")
	
	
	f1 = rt.TF1("f1",fitMassSpectrum,fmin,fmax,2)
	f1.SetParameters(10,-1.)
	hEBEB.Fit(f1, "0")
	
	
	#print "cumulative:"
	#hc = CumulativePlot(hEBEB, "down")
	#ctest = TCanvas("ctest","ctest",800,800)
	#ctest.SetLogx()
	#ctest.SetLogy()
	#ctest.cd()
	#hc.Draw("hist")
	#hEBEB.Draw("same ep")
	#print "// cumulative"
	
	# '''
	
	
	'''
	print "__EBEB"
	for binwidth in range(10, 41, 2):
		fmax = maxValFromBinwidth(fmin,maxvalue,binwidth)
		print fmax
		binWidthPlot(hEBEB,"EBEB",binwidth,fmin,fmax)
		raw_input()
	
	# '''
	
	# EBEE
	# 
	fmin = 330
	
	fmax = maxValFromBinwidth(fmin,maxvalue,binEB)
	
	hEBEE = rebin(hEBEE, range(fmin,fmax+1,binEE))
	hEBEE.GetXaxis().SetRangeUser(fmin, fmax)
	
	hEBEEc = CumulativePlot(hEBEE, "down")
	
	f2 = rt.TF1("f2",fitMassSpectrum,fmin,fmax,2)
	f2.SetParameters(10.,-1.)
	hEBEE.Fit(f2, "0")
	# '''
	
	'''
	print "__EBEE"
	for binwidth in range(10, 41, 2):
		fmax = maxValFromBinwidth(fmin,maxvalue,binwidth)
		print fmax
		#binWidthPlot(hEBEE,"EBEE",binwidth,fmin,fmax)
	
	# '''
	
	# cumulative plots for the fit model
	hf1 = CumulativePlot(createHistoFromFunction(f1, hEBEBc), "down")
	hf2 = CumulativePlot(createHistoFromFunction(f2, hEBEEc), "down")
	
	#create canvases with ratio plots
	rEBEB = canvasCreator(hEBEB, f1, hEBEB.GetName()+"_bin"+str(binEB))
	rEBEE = canvasCreator(hEBEE, f2, hEBEE.GetName()+"_bin"+str(binEE))
	
	rEBEB.draw()
	rEBEB.addHisto(hEBEBc, hf1)
	
	rEBEE.draw()
	rEBEE.addHisto(hEBEEc, hf2)
	# '''
	
	

def binWidthPlot(h,name,b,fmin,fmax):
	#
	htemp = rebin(h, range(fmin,fmax+1,b))
	htemp.GetXaxis().SetRangeUser(fmin, fmax)
	f = rt.TF1("f_"+name+str(b),fitMassSpectrum,fmin,fmax,2)
	f.SetParameters(10,-1.)
	htemp.Fit(f, "0")
	r = canvasCreator(htemp, f, name+"_bin"+str(b))
	r.draw()
	return 0


def maxValFromBinwidth(start,maximum,binwidth):
	# returns the maximum value where to go
	n = float(maximum-start)/float(binwidth)
	n = int(n)
	return start+binwidth*n
	


def canvasCreator(h, f, name):
	r = ratioCanvasHF(h,f,name)
	return r


# //////////////////////////////////////////////////////////////////////////////////////////////////
def main():
	_start = datetime.datetime.now()
	
	# write timestamp file
	if(FilelistDatetime()):
		print("file_dates.txt successfully written")
	else:
		print("file_dates.txt failed")
	
	rt.gROOT.SetBatch(rt.kTRUE)	# dont show the canvases
	
	
	print sys.argv
	
	
	# 
	massSpectrum()
	
	
	




if __name__ == "__main__":
	main()
	print "__main__ ... done!"
