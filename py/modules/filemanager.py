# 


# file and data management methods

import numpy as np
import datetime
import calendar
import sys
import os
from array import array








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

