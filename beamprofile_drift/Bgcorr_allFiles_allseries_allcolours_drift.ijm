run("Bio-Formats Macro Extensions");

//the raw (uncorrected files) need to be in a separte folder
rawFolder = getDirectory("Choose directory containing raw files");
//and will be saved in a separate folder that needs to exist
outputFolder = getDirectory("Choose a directory to save corrected files");
//read all files in Folder
files = getFileList(rawFolder);





//@ int darkframe = 2348;
#@ Integer (label="Darkframe value:",value=2348) darkframe



run("Set Measurements...", "min redirect=None decimal=0");
for (i = 0; i < files.length;i++) 
{
//reset all properties for every file   
interleaved = 0;
group = 0;
channels = 0;
frames = 0;
z = 0;
series = 0;
Ext.openImagePlus(rawFolder+files[i]);
Ext.setId(rawFolder+files[i]);
IJ.log(files[i]);
IJ.log(series);
Ext.getSeriesCount(series)
Ext.isInterleaved(interleaved)
Ext.isGroupFiles(group)
Ext.getSizeC(channels)
Ext.getSizeT(frames)
Ext.getSizeZ(z)
//IJ.log('format='+nseries);
//IJ.log('SeiresName = '+SeriesName);
//IJ.log('Interleaved? '+interleaved);
//IJ.log('grouped? '+group);
IJ.log('number of channels: '+channels);
IJ.log('number of z dimensions: '+z);
//IJ.log('number of frames: '+frames);
//IJ.log('z='+z)
IJ.log('seriesCount=' +series)
nseries = series;

//if more than 6 colors are corrected more colors need to added here...
colors = newArray('Green','Magenta','Yellow','Red','Cyan','Blue');

for (s = 1; s < nseries+1;s++)
{	
	series="series_"+s;
IJ.log("series =" + series);
IJ.log("Open " + files[i] + series);
	//display metadata and save it:  display_metadata!!
	run("Bio-Formats", "open=" + rawFolder+files[i] + " color_mode=Grayscale view=Hyperstack stack_order=XYCZT "+series );
	//waitForUser("see for data");
	rename("temp1");
	if (channels>1) {				//colors are only corrected if it exists 
		run("Split Channels");
		//selectWindow("C3-temp1");
		//close();
		
		
		 
		windows = getList("image.titles");
		IJ.log('There are'+windows.length+'images open')
		//windows = WindowManager.getImageTitles();
	
	
	// the first image is the original, the following ones the color channels. loop through them beginnning with the second.
	for(c = 1; c < (windows.length);c++){
			IJ.log(c)
			IJ.log(windows[c]);
			selectWindow(windows[c]);
			IJ.log('correct'+windows[c]);
			run("32-bit");
			run("Z Project...", "projection=[Average Intensity]");
			run("Gaussian Blur...", "sigma=80");
			run("Subtract...", "value=2348");
			run("Measure");
			max =(getResult("Max",0));
			run("Divide...", "value=max");
			selectWindow(windows[c]);
			run("Subtract...", "value="+darkframe+" stack");
			imageCalculator("Divide 32-bit stack", windows[c],"AVG_"+windows[c]);
			selectWindow("AVG_"+windows[c]);
			close();			
			selectWindow(windows[c]);
			setMinAndMax(0, 65535);
			run("16-bit");
			run(colors[c-1]);
			selectWindow("Results");
			run("Close");
		}
	//now the drift correction for every color:
	IJ.log("driftcorr");
	//Driftcorrection
	for(c = 1; c < (windows.length);c++){
		selectWindow(windows[0]);
		originalID = getImageID();
		setBatchMode(true);

		getSelectionBounds(sx, sy, sw, sh);
		//loop through slices, start at 2nd slice
		for (j = 2; j <= nSlices; j++) {
			selectImage(originalID);
			//go back to determine maximum at slice 1
			setSlice(1);
			makeRectangle(sx, sy, sw, sh);
			run("FFT");
			rename("a");
			selectImage(originalID);
			setSlice(j);
			makeRectangle(sx, sy, sw, sh);
			run("FFT");
			rename("b");
			run("FD Math...", "image1=a operation=Correlate image2=b result=c do");
			List.setMeasurements;
			cx = getWidth / 2;
			cy = getHeight / 2;
			max = 0;
			// the maximum should be somewhere in the center
			for (y = cy - 100; y <= cy + 100; y++) {
				for (x = cx - 100; x <= cx + 100; x++) {
					pixel = getPixel(x, y);			
					if (pixel > max) {
						max = pixel;
						dx = x;
						dy = y;
						}
					}
				}
			dx -= cx;
			dy -= cy;
			setResult("dx", j-2, dx);
			setResult("dy", j-2, dy);
	
			// close all temporary images
			selectImage("a");
			close();
			selectImage("b");
			close();
			selectImage("c");
			close();	
			}
			
	selectWindow(windows[c]);
	setBatchMode(false);
	originalID = getImageID();
	setBatchMode(true);
	getSelectionBounds(sx, sy, sw, sh);
	for (j = 2; j <= nSlices; j++) {
		selectImage(originalID);
		setSlice(j);
		run("Select All");
		dx=getResult("dx", j-2);
		dy=getResult("dy", j-2);
		//now translate the slices
		run("Translate...", "x=" + dx + " y=" + dy + " interpolation=Bilinear slice");
		setBatchMode(false);

	}
	}



//selectWindow("C2-temp1");
//	setBatchMode(false);
//
	//original = getImageID();
	//setBatchMode(true);

//	getSelectionBounds(sx, sy, sw, sh);

//	for (j = 2; j <= nSlices; j++) {
//		selectImage(original);
//		setSlice(j);
//		run("Select All");
//		dx=getResult("dx", j-2);
//		dy=getResult("dy", j-2);
//		run("Translate...", "x=" + dx + " y=" + dy + " interpolation=Bilinear slice");
//	}

//setBatchMode(false);
IJ.log("driftcorr finished");
	run("Merge Channels...", "c1=C1-temp1 c2=C2-temp1 create");
	
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	filename=("Corrected_"+files[i]);
	IJ.log("save : " + outputFolder+ "corrected_"+series+files[i]);

	saveAs("Tiff", outputFolder +  "corrected_"+files[i]+"_"+series);

	run("Close All");
	}
	else {


IJ.log("bgcorr");
	selectWindow("temp1");
	run("32-bit");
	run("Z Project...", "projection=[Average Intensity]");
	run("Gaussian Blur...", "sigma=80");
	run("Subtract...", "value=325");
	run("Measure");
IJ.log(i);
		max =(getResult("Max",0));
IJ.log("max =" + max);
		run("Divide...", "value=max");

	selectWindow("temp1");
	run("Subtract...", "value=340 stack");
	imageCalculator("Divide 32-bit stack", "temp1","AVG_temp1");
	selectWindow("AVG_temp1");
	close();
	selectWindow("temp1");
	setMinAndMax(0, 65535);
	run("16-bit");
	selectWindow("Results");
	run("Close");
	
	


IJ.log("driftcorr");
	//Driftcorrection
	selectWindow("temp1");
	original = getImageID();
	setBatchMode(true);

	getSelectionBounds(sx, sy, sw, sh);

	for (j = 2; j <= nSlices; j++) {
		selectImage(original);
		setSlice(1);
		makeRectangle(sx, sy, sw, sh);
		run("FFT");
		rename("a");
	
		selectImage(original);
		setSlice(j);
		makeRectangle(sx, sy, sw, sh);
		run("FFT");
		rename("b");
	
		run("FD Math...", "image1=a operation=Correlate image2=b result=c do");
	
		List.setMeasurements;
	
		cx = getWidth / 2;
		cy = getHeight / 2;
		max = 0;
	
		// the maximum should be somewhere in the center
		for (y = cy - 100; y <= cy + 100; y++) {
			for (x = cx - 100; x <= cx + 100; x++) {
				pixel = getPixel(x, y);
			
				if (pixel > max) {
					max = pixel;
					dx = x;
					dy = y;
				}
			
			}
		}
	
		dx -= cx;
		dy -= cy;
		setResult("dx", j-2, dx);
		setResult("dy", j-2, dy);
	
		// close all temporary images
		selectImage("a");
		close();
	
		selectImage("b");
		close();
	
		selectImage("c");
		close();	
	}
	selectWindow("temp1");
	setBatchMode(false);

	original = getImageID();
	setBatchMode(true);

	getSelectionBounds(sx, sy, sw, sh);

	for (j = 2; j <= nSlices; j++) {
		selectImage(original);
		setSlice(j);
		run("Select All");
		dx=getResult("dx", j-2);
		dy=getResult("dy", j-2);
		run("Translate...", "x=" + dx + " y=" + dy + " interpolation=Bilinear slice");
	}

setBatchMode(false);


IJ.log("driftcorr finished");

	
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	filename=("Corrected_"+files[i]);
	IJ.log("save : " + outputFolder+ "corrected_"+series+files[i]);

	saveAs("Tiff", outputFolder +  "corrected_"+files[i]+"_"+series);

	run("Close All");
	
	
	}
	
}
}

IJ.log("Finished!!");
