
import sys
sys.path.insert(0,"../../py")

import os
import datetime
import calendar
import time

import ROOT as rt
import numpy as np


import subprocess as sp

import inspect as ins

#from modules import lib_plot as lp

# ptest.py

def mystack():
	
	
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


class Foo(object):
	# test class
	def __init__(self):
		self.bar = None
	
	def __enter__(self):
		if self.bar != 'open':
			print 'opening the bar'
			self.bar = 'open'
		return self # this is bound to the `as` part
	
	def close(self):
		if self.bar != 'close':
			print 'closing the bar'
			self.bar = 'close'
	
	def __exit__(self, *err):
		self.close()
	
	
	def printme(self):
		print "hallo welt"


def fooCall():
		
		with Foo() as foo:
			a = foo
		
		a.printme()
		print a.bar


def multipardaughter(x, *y):
	
	if "hallo" in y:
		print "ja"
	
	if "dies" in y:
		print "auch"
		
	if y == []:
		print "leer"
	
	if y:
		print "jo!"
	else:
		print "no!"
	


def multipar():
	
	#multipardaughter(3,"a","hallo","dies","das",5)
	
	multipardaughter(3)
	
	


def replacement():
	
	s = "halloPTundMETundHToderTREFF"
	
	ss = parseTitle(s)
	
	print ss
	
	

def parseTitle(string):
    """
    replaces some keywords in string
    """
    replacements={
        'PT':'p_{T} [GeV]'
        ,'MET':'E_{T}^{miss} [GeV]'
        ,'HT':'H_{T} [GeV]'
        ,'HTG':'H_{T}^{gen} [GeV]'
        ,'MET':'E_{T}^{miss} [GeV]'
        ,'TREFF':'Trigger Effciency'
    }
    return reduce(lambda x, y: x.replace(y, replacements[y]), replacements,string)


def rangetest(start,maximum,binwidth):
	# 
	
	start = 320
	m = 2100
	
	b = 22
	
	n = (m-start)/b # int value, floatet values cutted
	
	return start+b*n
	


def functionstuff():
	
	f = TF1("f","x+3*x*x",10.,40.)
	
	print f.GetName()
	


def histofromf(f):
	
	xmin = 1
	
	return xmin




def lists():
	
	l = {}
	
	for i in range(0,10):
		l[i] = i*i
	
	it = l.itervalues()
	
	print it
	


def createGenerator():
	
	mylist = range(3)
	for i in mylist:
		yield i*i
	


def generatorsdaughter():
	pass



def generators():
	# iterators, generators and stuff
	# 'yield'
	
	mylist = [x*x for x in range(4)]
	print "mylist iterator: "
	for i in mylist:
		print i
	
	
	mygenerator = (x*x for x in range(4))
	print "mygenerator: "
	for i in mygenerator:
		print i
	
	print "mygenerator try second time: "
	for i in mygenerator:
		print i
	
	print "only one time possible!"
	print mygenerator
	
	print "now try yield: "
	mygenerator = createGenerator() # create a generator
	print mygenerator # mygenerator is an object!
	
	for i in mygenerator:
		print i
		
	print "second time: "
	for i in mygenerator:
		print i
	
	
	



def frametest():
	
	cf = ins.currentframe()
	framefile = ins.getfile(cf)
	filename = os.path.basename(os.path.abspath(framefile))
	path = os.path.dirname(os.path.abspath(framefile))
	absfile = path+'/'+filename
	
	print "currentframe: ",cf
	print "filename: ",filename
	print "path: ",path
	print "absolut file: ",absfile
	





def parseaxis():
	
	"""
	parses the axis titles of both axes:
	* quantities with units
	* y-axis: BIN->' / <binwidth>'
	"""
	
	# parse qunatities with units
	xt='p_{T} [GeV]'
	yt='H_{T}^{gen} [GeV] BIN'
	
	
	# parse y-axis
	if "BIN" in yt:
		# find bin-width
		binwidth="%g"%5
		# binwidth="%g"%xax.GetBinWidth(1)
	
		# check if x-axis has a unit
		import re
		pattern = re.compile('\[.*\]')
		unit = pattern.findall(xt)
		if unit:
			assert(len(unit)==1) # not two units
			unit=" "+unit[0][1:-1] # extract unit between [...]
		else:
			unit=' units' # dimensionless
		binwidth+=unit
		yt=yt.replace("BIN"," / %s"%binwidth)
	
	print pattern
	print xt
	print yt



class mydict(dict):
	def add(self,key,val):
		self[key] = val

def dictionaries():
	
	k = None
	
	if k is None:
		print "k!"
	if k is not None:
		print "no k!"
	
	
	a = mydict()
	a.add("first", "hallo")
	print a


def systemtest():
	
	import platform
	
	print platform.dist()
	print platform.linux_distribution()
	
	print os.name
	
	print sys.platform
	
	print platform.system()
	
	print platform.release()
	


def main():
	
	print sys.argv
	
	
	#functionstuff()
	#fooCall()
	#replacement();
	
	
	#lists()
	#generators()
	
	#frametest()
	
	#parseaxis()
	
	
	#dictionaries()
	
	systemtest()
	
	
	
	return 0



if __name__ == "__main__":
	main()
	print "__main__ ... done!"
