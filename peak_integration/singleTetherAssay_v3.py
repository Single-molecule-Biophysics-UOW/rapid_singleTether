from ij import IJ, ImagePlus
from ij.plugin.frame import RoiManager
from ij.plugin.filter import Analyzer
from ij.plugin import HyperStackConverter
from IOFunctions import pyIO
from bgCorrClass import Bgcorr
from driftCorrClass import DriftCorr
from ij import WindowManager

import coloc
import time
import os
from ij.measure import ResultsTable
from java.lang.System import getProperty # should crash your IDE
import java.awt.Robot as Robot
import java.awt.event.KeyEvent as KeyEvent

developmentMode = True	#set this to false if you're just using the script!

def deleteSource():
	"""deletes the compiled source file in the install dir, to make sure it's recompiled after changes are made.
	   tool for development! Should otherwise not be called!
"""

	fijiDir=os.path.abspath(getProperty('user.dir')+"/")
	LibPath = fijiDir + '/jars/Lib/'
	print LibPath
	os.remove(LibPath + "IOFunctions$py.class")
	os.remove(LibPath + "bgCorrClass$py.class")

#initialize imageFile
start = time.time()
#File = IJ.getFilePath("Choose File")

inputFolder = IJ.getDirectory("Choose input directory")
OutputFolder = IJ.getDirectory("Choose a directory to save output")

files = [x for x in os.listdir(inputFolder) if x.endswith('.tif')==True]

#identify reaction and prereaction:
preReaction =[]
reaction = []
for f in files:
	if 'prereaction' in f:
		preReaction.append(f)
	else:
		reaction.append(f)

print(preReaction,reaction)
for preFile,reactFile in zip(preReaction,reaction):	
	#File = inputFolder+filename
	#outname = (OutputFolder+filename[:-4]+'_corr.tif')
	#seriesCount,colorCount = pyIO.initializeFile(File)
	
	preIm = pyIO.loadFile(inputFolder+preFile,0)[0]
	reactIm = pyIO.loadFile(inputFolder+reactFile,0)[0]
	preFrames = preIm.getDimensions()[4]
	preChannels = preIm.getDimensions()[2]

	reactFrames = reactIm.getDimensions()[4]
	reactChannels = reactIm.getDimensions()[2]
	#corrIm = im.duplicate()
	
	IJ.log("image for both reaction and prereaction loaded,"+str(time.time()-start)+'seconds')


	#now find DNA peaks in reaction image and integrate in both channels
	#image= IJ.getImage()
	reactChannels = coloc.imSplitter(reactIm)
	if RoiManager.getInstance2() == None:
		RM = RoiManager()
	else:
		RM = RoiManager.getInstance2()
		RM.reset()
	reactIm.close()
	
	title = reactChannels[0].getTitle()
	#run maximum z-projection for finding all molecules
	IJ.run(reactChannels[0],"Z Project...", "stop=10 projection=[Max Intensity]")
	maxTitle = "MAX_"+title
	IJ.selectWindow(maxTitle)
	#run peak finder
	IJ.run("Peak Finder", "use_discoidal_averaging_filter threshold=2 selection_radius=3 minimum_distance=8")
	#save rois, the folder has to exist!
	RM.runCommand("Save", OutputFolder+reactFile[:-3]+'.zip')
	#max projection is not needed anymore!
	maxproj = WindowManager.getWindow(maxTitle)
	maxproj.close()
	IJ.run("Set Measurements...", "centroid redirect=None decimal=0")
	titles = [reactChannels[0].getTitle(),reactChannels[1].getTitle()]
	for i in range(len(reactChannels)):		#loop over all color channels and integrate trajectories
		IJ.run(reactChannels[i],"Integrate Trajectory", "background=2 integrate=[all slices] toImageJResultsTable")
		tableWindow = WindowManager.getWindow("Trajectory Integrations")	
		#IJ.run("Sort", "table=[Trajectory Integrations] column=slice group=trajectory ascending");
		tableWindow.setTitle(reactChannels[i].getTitle())	
		reactChannels[i].close()
	
	IJ.log("trajectories integrated in a total of {} seconds.".format(time.time()-start))

	#now integrate the same rois in the preReaction recording:
	preChannels = coloc.imSplitter(preIm)
	#original not needed anymore
	preIm.close()
	title = preChannels[0].getTitle()	
	titles = [title.getTitle() for title in preChannels]
	#preChannels[0].getTitle(),preChannels[1].getTitle()]
	for i in range(len(preChannels)):		#loop over all color channels and integrate trajectories
		IJ.run(preChannels[i],"Integrate Trajectory", "background=2 integrate=[all slices] toImageJResultsTable")
		tableWindow = WindowManager.getWindow("Trajectory Integrations")	
		#IJ.run("Sort", "table=[Trajectory Integrations] column=slice group=trajectory ascending");
		tableWindow.setTitle(preChannels[i].getTitle())	
		preChannels[i].close()		

if developmentMode == True:
	deleteSource()