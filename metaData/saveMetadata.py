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


inputFolder = IJ.getDirectory("Choose input directory")
OutputFolder = IJ.getDirectory("Choose a directory to save output")

files = [x for x in os.listdir(inputFolder) if x.endswith('.nd2')==True]


for filename in files:
	File = inputFolder+filename
	outname = (OutputFolder+filename[:-4])
	seriesCount,colorCount = pyIO.initializeFile(File)		
	IJ.run("Bio-Formats", "open="+File+" color_mode=Grayscale display_metadata rois_import=[ROI manager] view=[Metadata only] stack_order=Default");
	metadata = WindowManager.getActiveTable()
	IJ.saveAs("Text",outname+"_meta.csv")
	metadata.close()