from ij import IJ, ImagePlus
from ij.plugin.frame import RoiManager
from ij.plugin.filter import Analyzer
from ij.plugin import HyperStackConverter
from IOFunctions import pyIO

import table.ResultsTableUtil as tableUtil
from ij import WindowManager

import coloc
import time
import os
from ij.measure import ResultsTable
from java.lang.System import getProperty # should crash your IDE
import java.awt.Robot as Robot
import java.awt.event.KeyEvent as KeyEvent

developmentMode = True	#set this to false if you're just using the script!

def saveResultsTable(path):
	allTitles = tableUtil.getResultsTableTitles()
	#print(allTitles)
	title = tableUtil.getResultsTableTitles()[0]
	table = tableUtil.getResultsTable(title)
	table.save(path)	
	tableUtil.removeResultsTable(title)
def deleteSource():
	""" deletes the compiled source file in the install dir, to make sure it's recompiled after changes are made.
	   tool for development! Should otherwise not be called!
	"""

	fijiDir=os.path.abspath(getProperty('user.dir')+"/")
	LibPath = fijiDir + '/jars/Lib/'
	print LibPath
	os.remove(LibPath + "IOFunctions$py.class")
	os.remove(LibPath + "bgCorrClass$py.class")

#initialize imageFile


inputFolder = IJ.getDirectory("Choose input directory")
OutputFolder = IJ.getDirectory("Choose a directory to save output")

files = [x for x in os.listdir(inputFolder) if x.endswith('.tif')==True]

#
for f in files:	
	#File = inputFolder+filename
	#outname = (OutputFolder+filename[:-4]+'_corr.tif')
	#seriesCount,colorCount = pyIO.initializeFile(File)
	
	Im = pyIO.loadFile(inputFolder+f,0)[0]

	Frames = Im.getDimensions()[4]
	#Channels = preIm.getDimensions()[2]

	
	if RoiManager.getInstance2() == None:
		RM = RoiManager()
	else:
		RM = RoiManager.getInstance2()
		RM.reset()
	#reactIm.close()
	
	title = Im.getTitle()
	#run maximum z-projection for finding all molecules
	IJ.run(Im,"Z Project...", "stop=10 projection=[Max Intensity]")
	maxTitle = "MAX_"+title
	IJ.selectWindow(maxTitle)
	#run peak finder
	IJ.run("Peak Finder", "use_discoidal_averaging_filter threshold=2 selection_radius=3 minimum_distance=8")
	#save rois, the folder has to exist!
	#RM.runCommand("Save", OutputFolder+reactFile[:-3]+'.zip')
	#max projection is not needed anymore!
	maxproj = WindowManager.getWindow(maxTitle)
	maxproj.close()
	#IJ.run("Set Measurements...", "centroid redirect=None decimal=0")
	#titles = [reactChannels[0].getTitle(),reactChannels[1].getTitle()]
	#for i in range(len(reactChannels)):		#loop over all color channels and integrate trajectories
	#integrate trajectories
	
	IJ.run(Im,"Integrate Trajectory", "background=3 integrate=[all slices] toImageJResultsTable")
	tableWindow = WindowManager.getWindow("Trajectory Integrations")	
		#IJ.run("Sort", "table=[Trajectory Integrations] column=slice group=trajectory ascending");
	tableWindow.setTitle(Im.getTitle())	
	Im.close()
	saveResultsTable(OutputFolder+f[:-3]+'.csv')
	
	

	

if developmentMode == True:
	deleteSource()