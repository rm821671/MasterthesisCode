#

from ROOT import *
from lib_plot import *
from math import *
from array import array


def main():
	
	path = "/user/rmeyer/mergedNTuples/";
	ndata = 6
	datasets = [	"DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",
					"TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_v01.root",
					"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",
					"ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v01.root",
					"DoubleEG_Run2015D-05Oct2015-v1.root",
					"DoubleEG_Run2015D-PromptReco-v4.root"
				]
				
	
	

	
	#print 14.5e-15
	
	
	#print giveDataPath("DYJet")
	
	
	
	
	
	
	t = {"k": 3,
		"p": 5,
		"c": 2
		}
		
	
	a = dict.fromkeys(t.keys())
	
	
	
	print a
	
	print min(t, key=t.get) # returns the key corresponding to the minimum value
	
	
	print("try his = TH1F():")
	his = TH1F()
	his.SetName("hallo")
	print his.GetName()
	
	
	return True
	
	
	




if __name__ == "__main__":
	main()



