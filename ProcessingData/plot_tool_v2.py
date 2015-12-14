
from ROOT import *
from lib_plot import *
from math import *
from array import array



# global objects
csize = 800


def bw(x, pars):
	# breit wigner
	return pars[0]*1./(pow((pars[1]-x[0]),2) + pow((x[0]*pars[2]),2))
	

def fit_tool(h, f, c):
	h.Scale(1./float(h.GetEntries()))
	h.Fit(f)	
	c.cd()
	f.Draw("same")
	
	return True


def fake_dep(fr, var):
	pass

def rawplot(datatag, lines=6, lineoff = 0):
	'''
	parameter:
		lines, gives the number of files written per run
		lineoff, means the offset; by default zero, so it starts with the "latest" lines
	'''
	rebin = 2
	
	datakeys = {	"DYJetsToLL": 0,
					"TTGJets": 2,
					"WGToLNuG": 4,
					"ZGTo2LG": 6,
					"Run2015D-05Oct2015": 0,
					"Run2015D-PromptReco": 1
				}
	defcolor = {	0: kGreen}
	
	# get number of lines from file:
	total_lines = inputFileNumberOfLines()
	h = {} #histograms
	
	stacks = {}
	
	# read the files:
	for ll in range(total_lines-lineoff-lines+1,total_lines-lineoff+1):
		# filename
		strInputFile = inputFile(ll)
		for key in datakeys:
			if(key in strInputFile):
				hname = key
		#print hname
		h[hname] = {}
		try:
			data = TFile(strInputFile)
		except:
			return False		
		# loop all histograms and search for raw ones
		for key in data.GetListOfKeys():
			kgn = key.GetName()
			if datatag in kgn:
				#print type(data.Get(kgn))
				h[hname][kgn] = data.Get(kgn)
				h[hname][kgn].SetDirectory(0)
			
		
		
		stacks = dict.fromkeys(h[hname].keys()) # THStack for mc
	
	l = dict.fromkeys(stacks.keys()) # legends
	dh = dict.fromkeys(stacks.keys()) # data histogram(s)
	added_data = dict.fromkeys(stacks.keys())
		
	#initialize stacks with THStacks:
	for key in stacks.keys():
		stacks[key] = THStack()
		l[key] = TLegend(0.9,0.7,0.65,0.9)
		dh[key] = {}
		
	
	c = {} # canvases
	cstyle = ""
	
	# scale montecarlo samples and add to stack
	for key in h.keys():
		if "Run2015" not in key:
			sf = giveMClumiScalefactor(key)
			for his in h[key].keys():
				#print key, his
				h[key][his].Scale(sf)
				h[key][his].Rebin(rebin)
				h[key][his].SetFillColor(defcolor[0]+datakeys[key])
				h[key][his].SetLineColor(defcolor[0]+datakeys[key])
				l[his].AddEntry(h[key][his], key)
				stacks[his].Add(h[key][his])
	
	
		if "Run2015D-05Oct2015" in key:
			for his in h[key].keys():
				h[key][his].Rebin(rebin)
				dh[his][key] = h[key][his]
				#dh[his][key].SetFillColor(kBlack)
				dh[his][key].SetMarkerStyle(20) # full square
				dh[his][key].SetMarkerColor(kBlack)
				l[his].AddEntry(h[key][his], key)
				#print "key, his: ", key, his
				#print "his, key: ", his, key
	
	
		if "Run2015D-PromptReco" in key:
			for his in h[key].keys():
				h[key][his].Rebin(rebin)
				dh[his][key] = h[key][his]
				#dh[his][key].SetFillColor(kBlue)
				#dh[his][key].SetMarkerStyle(23) # full triangle up
				#dh[his][key].SetMarkerColor(kBlue)
				#l[his].AddEntry(h[key][his], key)
	
	# add reco and prompt reco
	for key in dh.keys():
		dh[key]["Run2015D-05Oct2015"].Add(dh[key]["Run2015D-PromptReco"])
		
	
	
	for key in stacks:
		print "key: ", key
		#print "key in stacks: ", key
		c[key] = TCanvas(key,key,csize,csize)
		dh[key]["Run2015D-05Oct2015"].Draw("ep")
		stacks[key].Draw("same hist")
		dh[key]["Run2015D-05Oct2015"].Draw("same ep")
		l[key].Draw()
	
	
	# test of ratio
	rh = Ratio(		dh["diphotonraw_photon_pt_test"]["Run2015D-05Oct2015"], 
					dh["diphotonraw_photon_pt_test"]["Run2015D-PromptReco"]
				)
	print "run bins: ", dh["diphotonraw_photon_pt_test"]["Run2015D-05Oct2015"].GetNbinsX()
	print "promptreco bins: ", dh["diphotonraw_photon_pt_test"]["Run2015D-PromptReco"].GetNbinsX()
	
	#ct = TCanvas("ct","ct",csize,csize)
	#ctt = TCanvas("ctt","ctt",csize,csize)
	
	#rh.hisdraw()
	#rh.hisdraw(ctt)
	#rh.getHis().Draw("ep")
	#print "drawed"
	
	myc = myCanvasRatio(dh["diphotonraw_photon_pt_test"]["Run2015D-05Oct2015"],
						rh,
						0)
						
	myc._create()
	
	'''
	cv = TCanvas("cv","cv",csize,csize)
	tp = TPad("tp","tp",0,0.3,1,1)
	tp.SetBottomMargin(0)
	tp.Draw()
	tp.cd()
	rh.getHis().Draw("ep")
	'''
	
	raw_input()



def main():
	_start = datetime.datetime.now()
	
	# write timestamp file
	if(FilelistDatetime()):
		print("file_dates.txt successfully written")
	else:
		print("file_dates.txt failed")
	
	#gROOT.SetBatch(kTRUE)	# dont show the canvases
	
	c = {}		# canvases
	#csize = 700
	
	hnumber = 0
	hnames = {}	# names of the histograms
	
	hnameplot = {}	# 
	h = {}		# histograms 
	h2 = {}		# 2d histograms
	h2p = {}	# projected histograms
	
	temph2 = 0
	
	#print sys.argv
	
	
	#rawplot("rawdata",6)
	
	rawplot("diphotonraw",6)
	
	#rawplot("diphotonZpeak",6)
	
	
	
	
	
	"""
	
	strInputFile = inputFile() # get latest created file
	
	print "input: ", strInputFile
	
	try:
		data = TFile(strInputFile)
		print "file successfully loaded"
	except:
		print "file load failed"
		sys.exit(0)
	
	for key in data.GetListOfKeys(): # loop all keys
		hnames[hnumber] = key.GetName()
		#print hnames[hnumber]
		c[hnumber] = TCanvas("c"+str(hnumber),hnames[hnumber],csize,csize)
		
		if "rawdata" in hnames[hnumber]:
			
		
		
		if "h2" in hnames[hnumber]:			
			print hnames[hnumber], hnumber
			h2[temph2] = data.Get(key.GetName())
			h2[temph2].Draw("colz")
			temph2 += 1
		
		hnumber += 1
	
	
	
	for k in range(0,temph2):
		
		
	
	
	myList = TList()
	for key in c.keys():
		myList.Add(c[key])
	
	strMyFile = GiveOutputString(strInputFile)
	print 'output: ',strMyFile
	
	myFile = TFile(strMyFile,"RECREATE")
	gFile.WriteObject(myList,"list")	
	
	end_ = datetime.datetime.now()
	print "runtime: ", end_ - _start
	
	#"""



	# sdf


if __name__ == "__main__":
	main()
	print "__main__ ... done!"
