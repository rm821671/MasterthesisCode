
import datetime
import calendar
import time

from ROOT import *
import numpy as np
import sys, os

import subprocess as sp


from lib_plot import *

def main():
	
	
	#path = "../plots"
	today = datetime.date.today()
	print today
	
	#print today.isocalendar()
	
	#path = "../plots/bille/"
	
	folder = GivePlotfolder()
	
	strFile = inputFile()
	strOutFile = GiveOutputString(strFile)
	
	print "GivePlotFolder(): ", folder
	print "inputFile(): ",strFile
	print "GiveOutputString(inputFile()):",strOutFile
	
	#
	
	
	
	"""
	commit=sp.check_output(["git","log","-1","--pretty=format:%h"])
	status=sp.check_output(["git","status"])
	modified="modified" in status
	
	print "%s%s"%(commit, ("*" if modified else ""))
	
	#b = datetime.datetime.utcnow()

	#c = calendar.timegm(b.timetuple())
	#c = calendar.timegm(datetime.datetime.utcnow().timetuple())

	#s = 'jo_' + str(c) + '.root'

	#print d


	#s = "../selector_rootFiles/my_selector_results_mcDY_1446802738.root"
	
	strFile = inputFile()
	hnames = []
	print "file to load: ", strFile
	
	try:
		data = TFile(strFile)
		print "file load successful"
	except:
		print "file load failed"
		sys.exit(0)
	
	for key in data.GetListOfKeys():
		hnames.append(key.GetName())
	
	
	print hnames
	"""
	
	
	"""
	s = "../plot_rootFiles/my_plot_results_v1_1447706232.root"
	ss = s.split("_")
	d = int(ss[len(ss)-1].split(".")[0])
	
	print d
	print datetime.datetime.fromtimestamp(d)
	
	FilelistDatetime()
	
	fsize = os.path.getsize(s)
	print fsize, type(fsize)
	
	
	print inputFile(15)
	print inputFile()
	"""
	
	
	"""
	try:
		data = TFile(strFile)
	except:
		print "failed"
		sys.exit(0)
	
	
	
	lists = []
	list_count = 0
	h = []
	
	for it in data.GetListOfKeys():
		#print it.GetName()
		lists.append(data.Get(it.GetName()))
		h.append(list_count)
		for l in lists[list_count]:
			#print l.GetName()
			h[list_count].append(l.Get(l.GetName()))



	td = {"1":"a","2":"b","3":"c","4":"d"};



	print td.itervalues()

	"""




	"""
	print s,"\n",out

	myList = TList()
	n1 = TNamed("name1","title1")
	n2 = TNamed("name2","title2")


	myList.Add(n1)
	myList.Add(n2)

	myFIle = TFile("pytest.root","UPDATE")

	gFile.WriteObject(myList,"list")


	"""

if __name__ == "__main__":
	main()
	print "__main__ ... done!"
