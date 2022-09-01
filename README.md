# rapid_singleThether
All code used for the data analysis in "Rapid single-molecule characterisation of enzymes involved in nucleic-acid metabolism"

## ImageJ script: integrate_peaks_singlecolor.py

This script runs in ImageJ/FIJI [see imagej.net](https://imagej.net/)
Prerequisites:
* s_m_b.jar needs to be copied in your ``fiji.app/plugins`` folder.
* coloc.py needs to be copied in your ``fiji.app/jars/Lib`` folder. create it if it does not exist
* IOFunctions.py needs to be copied in your ``fiji.app/jars/Lib`` folder.

Now just open the scipt in IMageJ and click run.
 
The scipt will prompt the user to use two directories. One containing the image data. The results will be saved to the second one. The script will produce tables with the integrated intensities for each detected peak. These tables will be saved in the choosen folder with the same filename as the image data.
 
 ## integrate peaks in multi channel files:
 
 The same prerequisites as above.
 Integrates peaks in any number of channels, however peaks will only be detected in the first channel.
 
 ## phi29_integrate_peaks
 
 For the phi29 data we recorded image data before the enzymatic reaction was initiated. We did this to subtract the recorded intensity prior to the reaction as backgorund. We need to integrate the same peaks in two colors over 2 separate files.
To do so all files need to be saved in hte same folder. The prereaction file and the actual reaction file are assumed to have the same name apart from a prefix ``prereaction``.
Example: prereaction file: ``prereaction_1uMATP.tif``
            reaction file: ``reaction_1uMATP.tif``
The script will then produce 4 tables with integrated peaks in both colors and for reaction and prereaction.

## Python scripts

All other files are pure python scripts. They depend on ``rst_tools_v2`` and take the integrated peaks tables as inputs to
* calculate rate histograms of multi-turnover reactions: (``rates_by_piecewise_linear_fits.py``)
* calculate a synchronized average trajectory (``synchronized_mean.py``)
* rates of unbinding (see UvrD section in paper; ``unbinding_rate_single_file.py``)

## other dependencies:

The python scripts use a number of common open source python packages, such as numpy,pandas, scipy, sci-kit learn.
 
