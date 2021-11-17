from ij import IJ, ImagePlus
from ij.plugin.frame import RoiManager



imp = IJ.getImage()
stack = imp.getImageStack()
dim = imp.getDimensions()

proc = stack.getProcessor(1000)

stats = proc.getStats()


meanI = []
with open('C:/Users/shm975/Documents/tempData/211101/uvrD/meanI.csv','w') as f:
	f.write('frame,meanI\n')
	for frame in range(1,dim[2]):
		stats= stack.getProcessor(frame).getStats()
		meanI=stats.mean
		f.write(str(frame)+','+str(meanI)+'\n')

	



#IJ.run("Peak Finder", "use_discoidal_averaging_filter inner_radius=1 outer_radius=3 threshold=6 threshold_value=0 selection_radius=4 minimum_distance=8 background=40 slice");


