from ij import IJ, ImagePlus
from ij import WindowManager
from loci.formats import ImageReader


import time
import os
from ij.measure import ResultsTable
from java.lang.System import getProperty # should crash your IDE
import java.awt.Robot as Robot
import java.awt.event.KeyEvent as KeyEvent

def initializeFile(path):
	"""method for getting the series count and number color channels without opening the file"""
	IR = ImageReader()
	IR.setId(path)
	SeriesCount = IR.getSeriesCount()
	colorCount = IR.getSizeC()
	return SeriesCount,colorCount


inputFolder = IJ.getDirectory("Choose input directory")
OutputFolder = IJ.getDirectory("Choose a directory to save output")

files = [x for x in os.listdir(inputFolder) if x.endswith('.nd2')==True]


for filename in files:
	File = inputFolder+filename
	outname = (OutputFolder+filename[:-4])
	seriesCount,colorCount = initializeFile(File)		
	IJ.run("Bio-Formats", "open="+File+" color_mode=Grayscale display_metadata rois_import=[ROI manager] view=[Metadata only] stack_order=Default");
	windows = WindowManager.getAllNonImageWindows()
	for i in windows:
		print(i.getTitle())
		if "Metadata" in i.getTitle():			
			metadata=i
	 #= WindowManager.getActiveTable()
	try:
		IJ.saveAs("Text",outname+"_meta.csv")
		metadata.close()
	except:
		print("could not find the window with metadata, something is wrong.")
	