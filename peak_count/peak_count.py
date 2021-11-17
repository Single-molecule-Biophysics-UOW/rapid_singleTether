from ij import IJ, ImagePlus
from ij.plugin.frame import RoiManager
import os
from IOFunctions import pyIO
from ij import WindowManager
import table.ResultsTableUtil as tableUtil
import ij.measure.ResultsTable as rt
from ij.gui import Roi
import coloc

def saveResultsTable(path):
	allTitles = tableUtil.getResultsTableTitles()
	#print(allTitles)
	title = tableUtil.getResultsTableTitles()[0]
	table = tableUtil.getResultsTable(title)
	table.save(path)	
	tableUtil.removeResultsTable(title)
	
	
inputFolder = IJ.getDirectory("Choose input directory")
OutputFolder = IJ.getDirectory("Choose a directory to save output")

files = [x for x in os.listdir(inputFolder) if x.endswith('.tif')==True]
n=1
for f in files:
	image = pyIO.loadFile(inputFolder+f,0)[0]
	#Channels = coloc.imSplitter(image)
	#print(Channels)
	DNA = image#Channels[0]
	#Channels[1].close()
	im = DNA
	#maybe:
	roi = Roi(148, 140, 264,264)
	im.setRoi(roi)
	#print(f)
	stack = im.getImageStack()
	dim = im.getDimensions()
	#print(dim)
	#now find mean and standard deviation on last frame:
	lastframe = dim[3]
	stats = stack.getProcessor(lastframe).getStats()
	meanI=stats.mean
	stdI = stats.stdDev

	print('stats of orig:Mean: {},Std: {}'.format(meanI,stdI))
	#get ImagePlus of lastFrame:
	lastframeImp = ImagePlus('last Frame' ,stack.getProcessor(lastframe))
	#show it

	
	IJ.run(lastframeImp,"Discoidal Averaging Filter", "inner_radius=1 outer_radius=4 slice");
	
	
	stats = lastframeImp.getProcessor().getStats()
	meanI=stats.mean
	stdI = stats.stdDev
	print('stats of discoidal avg:Mean: {},Std: {}'.format(meanI,stdI))
	print("threshold={}".format(meanI+1.5*stdI))
	lastframeImp.close()
	#print(meanI+2*stdI)
	#print('Mean:{}, Std:{}'.format(meanI,stdI))

	IJ.run(im,"Peak Finder", "use_discoidal_averaging_filter inner_radius=1 outer_radius=3 threshold_value={} selection_radius=4 minimum_distance=6 just_count_peaks background=3 stack".format(meanI+1.5*stdI));
	#print(meanI,stdI)
	#windows = WindowManager.getAllNonImageWindows()
	#peak_count = WindowManager.getWindow('Peak_Count')
	#JFrameToCsv(peak_count)
	#IJ.selectWindow('Peak_Count')
	#IJ.save(peak_count,OutputFolder+f)
	#for i in windows:
	#	print(i.title)
	#	if i.title == 'Peak_Count':
	saveResultsTable(OutputFolder+f[:-3]+'.csv')
	print('finished {} out of {}'.format(n,len(files)))
	im.close()
	n+=1
	
