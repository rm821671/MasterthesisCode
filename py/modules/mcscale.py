#




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
	returns lumi in pb
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
