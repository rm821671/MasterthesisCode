# ptest.py

import sys
from ROOT import *


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
	




def main():
	
	print sys.argv
	
	
	#functionstuff()
	#fooCall()
	#replacement();
	
	
	lists()
	
	
	
	
	return 0



if __name__ == "__main__":
	main()
