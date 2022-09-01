import os
from ij.gui import Roi
from ij.plugin.frame import RoiManager
from ij.io import FileSaver
from java.awt import Rectangle
from ij import IJ, ImagePlus
from ij import plugin
from ij import WindowManager
from itertools import product


#
#IJ.run("Peak Finder", "uyrukrfintegrate=true")

def imSplitter(image):
	splitter = plugin.ChannelSplitter()
	image = splitter.split(image)
	return image
def surrounding(A,r,border = (0,512,0,512)):
    """
    INPUT: A, tuple of X,Y coordinates
            r, radius of surrounding
            optional: bounding rectangle as X0,Xmax,Y0,Ymax,
                      coordinates outside this rectangle will be omitted
    Returns: List of (X,Y) coordinates within radius
    """
    dX = [i for i in range(int(A[0]-r),int(A[0]+r+1),1) if i >= border[0] and i < border[1]]
    dY = [i for i in range(int(A[1]-r),int(A[1]+r+1),1) if i >= border[2] and i < border[3]]
    XY = product(dX,dY)
    return XY
def addListToRM(rois,RMinstance,reset=True):
	RM = RMinstance
	if reset==True:
		RM.reset()
	for i in rois:
		index = RM.getCount()
		RM.add(i,index)
def xyListToRoi(xyList,r):
	"""
	creates a rectangular Roi of width=w and height = h around tuple(x,y,t) and adds them to the RoiManager (without resetting it)
	returns list of Rois
	"""
	color = 0
	RM=RoiManager().getInstance()
	roiList = []
	w=2*r
	h=w
	for j in xyList:	
		xyRoi = Roi(j[0]-r,j[1]-r,w,h)		
		roiList.append(xyRoi)
	return roiList

if __name__ in ['__builtin__', '__main__']:

	#get Current image:
	image = IJ.getImage()
	#splitt channels:
	channels = imSplitter(image)


	#channels[0].show()
	#IJ.setActiveImage(channel[0])
	#	RM = RoiManager.getInstance()
	if RoiManager.getInstance2() == None:
		RM = RoiManager()
	else:
		RM = RoiManager.getInstance2()
		RM.reset()
	IJ.run("Set Measurements...", "centroid redirect=None decimal=0");

	"""
	#the amount of rois fond is a bit weird and different from the manual way still!
	"""

	#create list to save all rois in
	RoiDicts = []
	for i in range(len(channels)):
		#find peaks in channels
		#print("find rois in channel" + str(i) +"now \n")
		IJ.run(channels[i],"Peak Finder", "use_discoidal_averaging_filter threshold=1 selection_radius=4 minimum_distance=7");
		RM.runCommand("Measure")
		table = WindowManager.getWindow("Results").getResultsTable()
		#for some reason X and Y are columns number 6 and 7!
		#create dictionary containing roi number as keys and coordinates as values (Round to integers!):
		#create dictionary, by looping rows			 X					   Y						row number
		ChannelRoiDict = dict(((int(round(table.getValue(6,row))),int(round(table.getValue(7,row)))),row) for row in range(table.size()))
		RoiDicts.append(ChannelRoiDict)
		table.reset()
		#reset RoiManager to avaoid confusion!
		RM.reset()

	#print("table size channel1: {} table size channel2: {}".format(len(RoiDicts[0]),len(RoiDicts[1])))
	#print(RoiDicts[0].keys())
	#print("now the second channel:")
	#print(RoiDicts[1].keys())
	#now look for colocalization, loop through channel1 dict, and check if there is a peak detected in dict2

	"""
	The trick: In the second dictionary we will create a dictionary with the coordinates as keys! So asking if a roi exists with specific coordinates becomes much much faster!
	Just make sure to catch the error thrown when a key doesn't exist!
	"""
	print("number of rois in channel1: {}".format(len(RoiDicts[0].keys())))
	#loop over coordinates (keys of dict 1):
	coloc = []
	for i in RoiDicts[0].keys():
		#calculate surounding of this peak:
		s = surrounding(i,2)
		#now look if any of the surounding pixels is a peak in the other color channel:
		for j in s:
			try:
				RoiDicts[1][j]
				coloc.append(j)
				break #skip to next surrounding if one member was found to colocalize to avoid counting multiple rois per surounding
			except KeyError:
				continue
	
	#print(coloc)

	colocRois = xyListToRoi(coloc,4)

	print("colocalized Rois: {}".format(len(coloc)))

	addListToRM(colocRois,RM)





