# 

import ROOT as rt
from lib_filemanagement import *
import numpy as np
import datetime
import calendar
import sys
import os




def checkRebinningConsistence( axis, newBinning ):
	# gives a warning if the bin widths are not consistent with the bin range
	# i.e. range 100-200, width 30
	oldBinning = []
	for i in range(axis.GetNbins()+1):
		oldBinning.append( axis.GetBinUpEdge(i) )
	for i in newBinning:
		if i not in oldBinning: print "New bin edge is not compatible with old binning", i


def rebin( h, binEdges ):
	# usage: h = rebin( hold, range(340,2001,20) )
	checkRebinningConsistence(h.GetXaxis(), binEdges)
	binEdgesArr = array( 'd', binEdges )
	hname = h.GetName() + "_rb"
	hnew = h.Rebin( len(binEdges)-1, hname, binEdgesArr )
	return hnew
