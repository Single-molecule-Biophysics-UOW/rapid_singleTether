






import os
import javax.swing as swing
from java.awt import BorderLayout, GridLayout,Dimension
import java
from ij.io import Opener

from loci.formats import ImageReader
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import time

class pyIO():
	@staticmethod
	def saveData(image,path,closeIm=True):
		FS = FileSaver(image)	
		FS.saveAsTiff(path+'.tif')
		if closeIm == True:
			image.close()
	@staticmethod
	def initializeFile(path):
		IR = ImageReader()
		IR.setId(path)
		SeriesCount = IR.getSeriesCount()
		colorCount = IR.getSizeC()
		return SeriesCount,colorCount
	@staticmethod
	def loadFile(path,series):			
		options = ImporterOptions()
		options.setSeriesOn(series,True)
		options.setVirtual(True)	#something goes wrong when it is set to true....
		options.setId(path)
		im = BF.openImagePlus(options)
		return im