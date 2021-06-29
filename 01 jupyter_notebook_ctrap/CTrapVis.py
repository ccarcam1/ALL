"""
This script is the first version of a GUI designed to dynamically look through
.h5 files from a Lumicks C-Trap instrument. Users are able to choose what plots
they want (Force vs. Distance, Force vs. Time, RGB image, or a combination of 
RGB and one of the force options) with options to customize the plot to add a 
Title, scale RGB image values, shift axis of all available axis. 

If both plots are dependent on time the x-axes are linked by the manual entries,
but the individual plots can be zoomed in/out of using the interactive toolbar.
The save name for both the image and the metadata files can be customized by the 
"File Name to Save" entry and the "Image Format" option. The "Draw Plot" button
pulls in all of the GUI information and draws/plots the correct plot back on 
the interface. The "Show RGB Histogram" button shows the distribution of red, 
green, and blue pixel intensities after any image manipulation. The "Quit"
button destroyed the tkinter GUI and quits the python execution. The "Export 
Force For Origin" button pulls a separate Tkinter window to allow the user
to export desired force data to a .csv file for further analysis in software
packages like Origin. The "Export Image For CTrapViewer" exports .png files
compatible with Dr. Iodo Hellar's Lab CTrapViewer software. The "Export Image 
For ImageJ" button exports scan image in such a way that the axis are shown
as well as a timestamp of the image for easy transformation into an ImageJ
montage.

Tested in Python 3.7.6
Modules Used: numpy, lumicks.pylake, matplotlib, tkinter, pandas, glob, os
If any of the modules are not avaliable run "pip install moduleName" in terminal
or "conda install moduleName"

Author: John Watters (with significant borrowing from Michael Wasserman's code)
"""

# import necessary modules
import numpy as np
from lumicks import pylake
import matplotlib
#matplotlib.use("TkAgg")
matplotlib.use("Agg")
import os
import glob
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import pandas as pd

# Using a class to describe the GUI interface and provide the different functions
class CTrapGUI():
    """
    Upon initialization, build the core GUI system of labels and buttons with default values
    the "change folder" button can then be used to change the file directory to the file of your choice.
    """
    def __init__(self, master):         
        self.__version__ = "1.0.0"
        self.__versionDate__ = "9/26/20"
        self.__cite__ = "Watters, J.W. (2020) C-Trap .h5 Visualization GUI. Retrieved from https://harbor.lumicks.com/"

        """
        Function that when called takes all of the force option/plotting options from the interface
        and extracts the correct force and time/distance values from the .h5 file
        Used a different function for FD curves to limit confusion
        """
        def extract_force_relevant(h5file,timestampsForIndexing=('',''),timestampsForScanIndexing=('',''),multiScanShading=0):
            amtToDS = float(entryDownSample.get())
            forceString = forceChannelPulldown.get()
            distanceVarString = whichDistanceValue.get()
            downsampleOpt = checkValueDownsampleOpt.get()
            stringForceChannel = 'Force ' + forceString
            
            #get sample rate to determine the amount to downsample
            sample_rate = h5file['Force HF']['Force 1x'].sample_rate
            
            
            if multiScanShading == 0:
                """
                The combobox logical gate is to determine Force vs Time or Force vs. Distance.
                The len > 3 logical gate is to see if the user wants the weighted average of the X and Y forces on the bead or a singular channel.
                If downsampleOpt = 1 then the downsampled data will be extracted instead of the HF data.
                For the FD option --> only low frequency force data is acquired.
                """
                if comboboxForFTvsFD.get() == "Force-Time":
                    #force time data collection
                    if len(forceString) < 3:
                        #one channel option
                        if downsampleOpt == 1:
                            yData = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            xData = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                        else:
                            yData = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            xData = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps)
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                    else:
                        #average of channels for one of the beads
                        if downsampleOpt == 1:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yData = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            xData = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                        else:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            yData = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            xData = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps
                            maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                            xMin = np.min(xData)
                            xMax = np.max(xData)
                            xData = np.interp(xData, (xMin,xMax), (0, maxTime))
                            #assert 316 == 1
                else: #fd option
                    if len(forceString) < 3:
                        #get force data for a specific channel
                        yData = h5file["Force LF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        #get distance data that corresponds to the distance
                        xData = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                    else:
                        #get distance data
                        xComponent = h5file["Force LF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yComponent = h5file["Force LF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yData = np.sqrt(np.square(xComponent) + np.square(yComponent))
                        #get force data that corresponds to the distance
                        xData = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data                
                
                if forceString == '2x':
                    yData = yData * -1
                return yData, xData
            
            #logical gate for multiple scan images - to highlight what force regime you are in
            if multiScanShading != 0:
                maxTime = (int(timestampsForIndexing[1]) - int(timestampsForIndexing[0])) / 1e9
                maxScanTime = (int(timestampsForScanIndexing[1]) - int(timestampsForScanIndexing[0])) / 1e9
                #same FT vs FD gate
                if comboboxForFTvsFD.get() == "Force-Time":
                    #force time data collection
                    if len(forceString) < 3:
                        #one channel option
                        if downsampleOpt == 1:
                            #generate full data
                            yDataFull = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            xDataFull = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))                            
                            
                            #index into the Force file for the partial data to highlight
                            scanY = h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            scanX = (h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps)
                            xMin = np.min(scanX)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                        else:
                            yDataFull = h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            xDataFull = (h5file["Force HF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps)
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))                            
                            
                            #index into the Force file for the partial data to highlight
                            scanY = h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                            scanX = (h5file["Force HF"][stringForceChannel][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].timestamps)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                    else:
                        #average of channels for one of the beads
                        if downsampleOpt == 1:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yDataFull = np.sqrt(np.square(xComponent) + np.square(yComponent))

                            xDataFull = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).timestamps
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))

                            #partial index
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            yComponent = h5file["Force HF"]["Force " + forceString[0]+"y"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS)).data
                            scanY = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            scanX = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].downsampled_by(int(sample_rate/amtToDS))
                            scanX = scanX.timestamps
                            xMin = np.min(scanX)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                        else:
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"].data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"].data
                            yDataFull = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            xDataFull = h5file["Force HF"]["Force 1x"].timestamps
                            xMin = np.min(xDataFull)
                            xMax = np.max(xDataFull)
                            xDataFull = np.interp(xDataFull, (xMin,xMax), (0, maxTime))
                            
                            #partial index
                            xComponent = h5file["Force HF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            yComponent = h5file["Force HF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                            scanY = np.sqrt(np.square(xComponent) + np.square(yComponent))
                            scanX = h5file["Force HF"]["Force 1x"][timestampsForIndexing[0] : timestampsForIndexing[1]].timestamps
                            xMin = np.min(scanX)
                            scanX = (scanX - int(timestampsForIndexing[0])) / 1e9
                            
                else: #fd option
                    if len(forceString) < 3:
                        #get force data for a specific channel
                        yDataFull = h5file["Force LF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        #get distance data that corresponds to the distance
                        xDataFull = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        
                        #get partial data
                        scanY = h5file["Force LF"][stringForceChannel][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        scanX = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                    else:
                        #get force data
                        xComponent = h5file["Force LF"]["Force " + forceString[0] +"x"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yComponent = h5file["Force LF"]["Force " + forceString[0] +"y"][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        yDataFull = np.sqrt(np.square(xComponent) + np.square(yComponent))
                        #get distance data that corresponds to the force
                        xDataFull = h5file["Distance"][distanceVarString][timestampsForIndexing[0] : timestampsForIndexing[1]].data
                        
                        #get partial data
                        xComponent = h5file["Force LF"]["Force " + forceString[0] +"x"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                        yComponent = h5file["Force LF"]["Force " + forceString[0] +"y"][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                        scanY = np.sqrt(np.square(xComponent) + np.square(yComponent))
                        scanX = h5file["Distance"][distanceVarString][timestampsForScanIndexing[0] : timestampsForScanIndexing[1]].data
                
                if forceString == '2x':
                    yDataFull = yDataFull * -1
                    try:
                        scanY = scanY * -1
                    except:
                        pass
                    
                return yDataFull, xDataFull, scanY, scanX
        
        
        #code to modify RGB values
        def modify_rgb_image(RGB_code):
            bightnessAddition = int(entryBrightness.get())
            
            grayscaleOption = grayscaleOpt.get()
            
            
            if bugFixButton.get():
                RGB_shape = RGB_code.shape
                RGB_code = np.reshape(RGB_code,(RGB_shape[1],RGB_shape[0],3))
                
            
            if grayscaleOption == "None":
                mod_RGB = RGB_code[:,:,:] * [float(entryRed.get()), float(entryGreen.get()), float(entryBlue.get())]
                if bightnessAddition != 0:
                    mod_RGB = (mod_RGB + bightnessAddition)
            else:
                if grayscaleOption == "Red":
                     mod_RGB = RGB_code[:,:,0] * float(entryRed.get())
                elif grayscaleOption == "Green":
                    mod_RGB = RGB_code[:,:,1] * float(entryGreen.get())
                else:
                    mod_RGB = RGB_code[:,:,2] * float(entryBlue.get())
                
            mod_RGB = mod_RGB.astype(int)
            return mod_RGB
        
        """
        Code to reset the boundary entry values and to remove previous labels describing the maximums
        - max values are set to '-' to help other code recognize when a new plot is being made and new maxes need to be defined
        - min values are set to 0 - if you are analyzing a Force-Distance plot then this zero will be used as a logical gate to set
        the distance minimum value to the distance minimum boundary
        """
        def reset_axis_items(diffFileOpt):
            if diffFileOpt != 0:
                    entryTimeMin.delete(0, "end")
                    entryYForceMin.delete(0, "end")
                    entryDistMin.delete(0, "end")
                    entryYRGBMin.delete(0, "end")
                    entryTimeMax.delete(0, "end")
                    entryDistMax.delete(0,"end")
                    entryYForceMax.delete(0, "end")
                    entryYRGBMax.delete(0, "end")
                    entryScanWidthMin.delete(0,"end")
                    entryScanWidthMax.delete(0,"end")
                    entryTimeMin.insert(0,'0')
                    entryYForceMin.insert(0,'0')
                    entryDistMin.insert(0,'0')
                    entryYRGBMin.insert(0,'0')
                    entryScanWidthMin.insert(0,"0")
                    entryTimeMax.insert(0,'-')
                    entryYForceMax.insert(0,'-')
                    entryDistMax.insert(0,'-')
                    entryYRGBMax.insert(0,'-')
                    entryScanWidthMax.insert(0,"-")
                    
                    try:
                        labelTimeMax.grid_forget()
                    except:
                        pass
                    try:
                        labelDistMax.grid_forget()
                    except:
                        pass
                    try:    
                        labelYRGBMax.grid_forget()
                    except:   
                        pass
                    try:
                        labelYForceMax.grid_forget()
                    except:
                        pass
                    try:
                        labelScanWidthMax.grid_forget()
                    except:
                        pass
            return
        
        """
        Function that will be linked to both the Save Image and Draw Plot feature.
        This function resets the channel data and writes the channel data (writeMetadat(file)) 
            if the new file/component type is different from the last one plotted.
        The typeComponent GUI input is then scanned and then output is sent to either
            extractAndPlotKymo, extractAndPlotScan (and its subset extractAndPlotStack),
            and extractAndPlotFD.
        """
        def generateFigure():
            #remove previous metadata and add metadata description from the .h5 file
            def writeMetadata(file):
                global metadataLabel
                
                #remove the label of previous metadata
                try:
                    metadataLabel.grid_forget()
                except:
                    pass
                
                metadataText = file.description
                metadataLabel = tk.Label(frameForMetadata,text=metadataText,justify=tk.LEFT)
                metadataLabel.grid(row=1,column=0,rowspan=1,columnspan=1,sticky='w')
                return
            
            
            #Function to extract data type that is specific to the Kymo data type
            def extractAndPlotKymo(h5file,resetBoundsOpt):                
                global labelTimeMax
                global labelDistMax
                global labelYRGBMax
                global labelYForceMax
                global labelScanWidthMax
                
                
                kymoString = typePulldown.get()
                splitKymoString = kymoString.split('-')
                kymoNumber = splitKymoString[-1]
                kymoPointer = h5file.kymos[kymoNumber]
                
                reset_axis_items(resetBoundsOpt)
                
                minTimeIndex = kymoPointer.timestamps[0,0].astype(np.int64)
                maxTimeIndex = np.max(kymoPointer.timestamps).astype(np.int64)
                
                if plottingOpt.get() == "Both":
                    fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True)
                    
                    RGB_unaltered = kymoPointer.rgb_image
                    
                    dx = kymoPointer.json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    dt = (kymoPointer.timestamps[0,1]-kymoPointer.timestamps[0,0]) / 1000000000 #scan time in s
                    RGB_altered = modify_rgb_image(RGB_unaltered)
                    
                    maxTime = len(kymoPointer.timestamps[0,:])
                    numberPixels = len(RGB_unaltered)
                    maxTrueTime = maxTime*dt
                    maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                    if entryTimeMax.get() == '-':
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(maxTrueTime,3)))
                        labelTimeMax = tk.Label(frameForAxis,text=str(round(maxTrueTime,3)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryYRGBMax.get() == '-':
                        entryYRGBMax.delete(0, "end")
                        entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                        labelYRGBMax = tk.Label(frameForAxis,text=str(round(maxTrueDist,2)))
                        labelYRGBMax.grid(row=3,column=3)
                    
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (maxTrueTime / maxTrueDist))
                    
                    if grayscaleOpt.get() == "None":
                        ax1.imshow(RGB_altered, **{**default_kwargs})
                    else:
                        ax1.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
                    
                    ax1.set_ylabel(u'Position(\u03bcm)')
                    ax1.set_xlabel('Time(s)')
                    ax1.set_xlim([float(entryTimeMin.get()),float(entryTimeMax.get())])
                    ax1.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                    
                    if comboboxForFTvsFD.get() == "Force-Distance":
                        forceData, distData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                        ax2.plot(distData,forceData)
                        
                        if entryDistMax.get() == '-':
                            entryDistMin.delete(0, "end")
                            entryDistMin.insert(0,str(round(min(distData),2)))
                            entryDistMax.delete(0, "end")
                            entryDistMax.insert(0,str(round(max(distData),2)))
                            labelDistMax = tk.Label(frameForAxis,text=str(round(max(distData),2)))
                            labelDistMax.grid(row=2,column=3)
                        
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax2.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                        ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax2.set_ylabel('Force(pN)')
                        ax2.set_xlabel(u'Distance(\u03bcm)')
                        
                        forceString = forceChannelPulldown.get()
                        
                        if len(forceChannelPulldown.get()) > 2:
                            bead_ID = forceString[0]
                            combinedString = 'Bead ' + bead_ID + ' vs. Distance'
                        else:
                            combinedString = 'Channel ' + forceString + ' vs. Distance'
                        ax2.set_title(combinedString)
                        
                    else:
                        forceData, timeData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                        ax2.plot(timeData,forceData)

                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                        ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax2.set_ylabel('Force(pN)')
                        ax2.set_xlabel('Time(s)')
                        
                        forceString = forceChannelPulldown.get()
                        
                        if len(forceChannelPulldown.get()) > 2:
                            bead_ID = forceString[0]
                            if checkValueDownsampleOpt.get() == 1:
                                combinedString = 'Bead ' + bead_ID + ' Downsampled to ' + entryDownSample.get() + ' Hz'
                            else:
                                combinedString = 'Bead ' + bead_ID
                        else:
                            if checkValueDownsampleOpt.get() == 1:
                                combinedString = forceString + ' Channel Downsampled to ' + entryDownSample.get() + ' Hz'
                            else:
                                combinedString = forceString + ' Channel'
                        ax2.set_title(combinedString)

                        
                elif plottingOpt.get() == "Force Only":
                    fig, ax = plt.subplots(constrained_layout=True)
                    if comboboxForFTvsFD.get() == "Force-Distance":
                        forceData, distData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeIndex,maxTimeIndex))
                        ax.plot(distData,forceData)
                        
                        if entryDistMax.get() == '-':
                            entryDistMin.delete(0, "end")
                            entryDistMin.insert(0,str(round(min(distData),2)))
                            entryDistMax.delete(0, "end")
                            entryDistMax.insert(0,str(round(max(distData),2)))
                            labelDistMax = tk.Label(frameForAxis,text=str(round(max(distData),2)))
                            labelDistMax.grid(row=2,column=3)
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                        ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax.set_ylabel('Force(pN)')
                        ax.set_xlabel(u'Distance(\u03bcm)')
                        
                        forceString = forceChannelPulldown.get()

                    else:
                        minTimeIndex = kymoPointer.timestamps[0,0].astype(np.int64)
                        maxTimeIndex = np.max(kymoPointer.timestamps).astype(np.int64)
                        forceData, timeData = extract_force_relevant(h5file,timestampsForIndexing= (minTimeIndex,maxTimeIndex))
                        ax.plot(timeData,forceData)
                        
                        if entryTimeMax.get() == '-':
                            entryTimeMax.delete(0, "end")
                            entryTimeMax.insert(0,str(round(max(timeData),3)))
                            labelTimeMax = tk.Label(frameForAxis,text=str(round(max(timeData),3)))
                            labelTimeMax.grid(row=1,column=3)
                        
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                        ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax.set_ylabel('Force(pN)')
                        ax.set_xlabel('Time(s)')
                        
                        forceString = forceChannelPulldown.get()
                else:
                    RGB_unaltered = kymoPointer.rgb_image
                    dx = kymoPointer.json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    dt = (kymoPointer.timestamps[0,1]-kymoPointer.timestamps[0,0]) / 1000000000 #scan time in s
                    RGB_altered = modify_rgb_image(RGB_unaltered)
                    maxTime = len(kymoPointer.timestamps[0,:])
                    numberPixels = len(RGB_unaltered)
                    maxTrueTime = maxTime*dt
                    maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                    if entryTimeMax.get() == '-':
                        entryTimeMax.delete(0, "end")
                        entryTimeMax.insert(0,str(round(maxTrueTime,3)))
                        labelTimeMax = tk.Label(frameForAxis,text=str(round(maxTrueTime,3)))
                        labelTimeMax.grid(row=1,column=3)
                    if entryYRGBMax.get() == '-':
                        entryYRGBMax.delete(0, "end")
                        entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                        labelYRGBMax = tk.Label(frameForAxis,text=str(round(maxTrueDist,2)))
                        labelYRGBMax.grid(row=3,column=3)
                    
                    fig, ax = plt.subplots()
                    default_kwargs = dict(extent=[0, maxTrueTime, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (maxTrueTime / maxTrueDist))
                    
                    if grayscaleOpt.get() == "None":
                        ax.imshow(RGB_altered, **{**default_kwargs})
                    else:
                        ax.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
                    
                    ax.set_ylabel(u'Position(\u03bcm)')
                    ax.set_xlabel('Time(s)')
                    ax.set_xlim([float(entryTimeMin.get()),float(entryTimeMax.get())])
                    ax.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                return fig
            
            
            """
            This function is an offshoot of the extractAndPlotScans if the RGB data is 3 dimensional (a stack)
            The multiScanShading option tells the program whether or not to highlight the force data from the different scan image (vs. the stack)
            Currently multiScanShading option has not been tested for Force-Distance curves
            """
            def extractAndPlotMultipleScans(h5file,scanPointer,resetBoundsOpt):
                global labelTimeMax
                global labelDistMax
                global labelYRGBMax
                global labelYForceMax
                global labelScanWidthMax
                
                scanString = typePulldown.get()
                splitScanString = scanString.split('-')
                scanNumber = splitScanString[-1]
                scanPointer = h5file.scans[scanNumber]
                
                #plot options are already reset in the extractAndPlotStack header
                #see if the stack being observed is a new stack object
                if resetBoundsOpt == 1:
                    scanNumber = scanPointer.num_frames
                    scaleForStack.config(from_=1, to=scanNumber)
                    scaleForStack.set(1)
                    
                #scanPointerFrame = scanPointer(frame=2)
                timestampArray = scanPointer.timestamps
                scanTimeStamp = (timestampArray[int(scaleForStack.get())-1,:,:]) 
                stackRGB = scanPointer.rgb_image[int(scaleForStack.get())-1,:,:,:]
                minTimeScanObj = np.min(timestampArray[0,0,0]).astype(np.int64)
                maxTimeScanObj = np.max(timestampArray[int(scanPointer.num_frames)-1,:,:]).astype(np.int64)
                minTimeScan = scanTimeStamp[0,0].astype(np.int64)
                maxTimeScan = np.max(scanTimeStamp[:,:]).astype(np.int64)
                
                if plottingOpt.get() == "Both":
                    fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True)
                    RGB_unaltered = stackRGB
                    dx = scanPointer.json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    totalScanWidth = scanPointer.json['scan volume']['scan axes'][0]['scan width (um)']
                    dt = (scanPointer.timestamps[0,0,1]-scanPointer.timestamps[0,0,0]) / 1000000000 #scan time in s
                    RGB_altered = modify_rgb_image(RGB_unaltered)
                    
                    #maxTime = len(timestampArray[0,:])
                    numberPixels = len(RGB_unaltered)
                    #maxTrueTime = maxTime*dt
                    maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                    if entryScanWidthMax.get() == '-':
                        entryScanWidthMax.delete(0,"end")
                        entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                        labelScanWidthMax = tk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                        labelScanWidthMax.grid(row=5,column=3)
                    
                    if entryYRGBMax.get() == '-':
                        entryYRGBMax.delete(0, "end")
                        entryYRGBMax.insert(0,str(round(maxTrueDist,3)))
                        labelYRGBMax = tk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                        labelYRGBMax.grid(row=3,column=3)
                    
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                    
                    ax1.imshow(RGB_altered, **{**default_kwargs})
                    ax1.set_ylabel(u'Position(\u03bcm)')
                    ax1.set_xlabel(u'Width(\u03bcm)')
                    ax1.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                    ax1.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                    
                    if comboboxForFTvsFD.get() == "Force-Distance":
                        if multiScanPlotOpt.get() == 1:
                            forceData, distData, scanForceData, scanDistData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj),timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                            ax2.plot(distData,forceData)
                            ax2.plot(scanDistData, scanForceData)
                        else:
                            forceData, distData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                            ax2.plot(distData,forceData)
                        
                        if entryDistMax.get() == '-':
                            entryDistMin.delete(0, "end")
                            entryDistMin.insert(0,str(round(min(distData),2)))
                            entryDistMax.delete(0, "end")
                            entryDistMax.insert(0,str(round(max(distData),2)))
                            labelDistMax = tk.Label(frameForAxis,text=str(round(max(distData),2)))
                            labelDistMax.grid(row=2,column=3)
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax2.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                        ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax2.set_ylabel('Force(pN)')
                        ax2.set_xlabel(u'Distance(\u03bcm)')
                        
                        forceString = forceChannelPulldown.get()
                        
                        if len(forceChannelPulldown.get()) > 2:
                            bead_ID = forceString[0]
                            combinedString = 'Bead ' + bead_ID + ' vs. Distance'
                        else:
                            combinedString = 'Channel ' + forceString + ' vs. Distance'
                        ax2.set_title(combinedString)
                        
                    else:
                        if multiScanPlotOpt.get() == 1:
                            forceData, timeData, scanForceData, scanTimeData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj), timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                        
                            ax2.plot(timeData,forceData)
                            ax2.plot(scanTimeData, scanForceData)
                        else:
                            forceData, timeData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                            ax2.plot(timeData,forceData)
                        
                        if entryTimeMax.get() == '-':
                            entryTimeMin.delete(0,"end")
                            entryTimeMin.insert(0,'0')
                            entryTimeMax.delete(0,"end")
                            entryTimeMax.insert(0,str(round(max(timeData),3)))
                            labelTimeMax = tk.Label(frameForAxis,text=str(round(max(timeData),3)))
                            labelTimeMax.grid(row=1,column=3)
                            
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0,"end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                        ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax2.set_ylabel('Force(pN)')
                        ax2.set_xlabel('Time(s)')
                        
                        forceString = forceChannelPulldown.get()
                        
                        if len(forceChannelPulldown.get()) > 2:
                            bead_ID = forceString[0]
                            if checkValueDownsampleOpt.get() == 1:
                                combinedString = 'Bead ' + bead_ID + ' Downsampled to ' + entryDownSample.get() + ' Hz'
                            else:
                                combinedString = 'Bead ' + bead_ID
                        else:
                            if checkValueDownsampleOpt.get() == 1:
                                combinedString = forceString + ' Channel Downsampled to ' + entryDownSample.get() + ' Hz'
                            else:
                                combinedString = forceString + ' Channel'
                        ax2.set_title(combinedString)

                        
                elif plottingOpt.get() == "Force Only":
                    fig, ax = plt.subplots(constrained_layout=True)
                    if comboboxForFTvsFD.get() == "Force-Distance":
                        if multiScanPlotOpt.get() == 1:
                            forceData, distData, scanForceData, scanDistData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj),timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                            ax.plot(distData,forceData)
                            ax.plot(scanDistData, scanForceData)
                        else:
                            forceData, distData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                            ax.plot(distData,forceData)
                        
                        if entryDistMax.get() == '-':
                            entryDistMin.delete(0, "end")
                            entryDistMin.insert(0,str(round(min(distData),2)))
                            entryDistMax.delete(0, "end")
                            entryDistMax.insert(0,str(round(max(distData),2)))
                            labelDistMax = tk.Label(frameForAxis,text=str(round(max(distData),2)))
                            labelDistMax.grid(row=2,column=3)
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                        ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax.set_ylabel('Force(pN)')
                        ax.set_xlabel(u'Distance(\u03bcm)')
                        
                    else:
                        if multiScanPlotOpt.get() == 1:
                            forceData, timeData, scanForceData, scanTimeData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj), timestampsForScanIndexing=(minTimeScan, maxTimeScan),multiScanShading=1)
                        
                            ax.plot(timeData,forceData)
                            ax.plot(scanTimeData, scanForceData)
                        else:
                            forceData, timeData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj,maxTimeScanObj))
                            ax.plot(timeData,forceData)
                        
                        if entryTimeMax.get() == '-':
                            entryTimeMax.delete(0, "end")
                            entryTimeMax.insert(0,str(round(max(timeData),3)))
                            labelTimeMax = tk.Label(frameForAxis,text=str(round(max(timeData),3)))
                            labelTimeMax.grid(row=1,column=3)
                        
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=4)
                        
                        ax.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                        ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax.set_ylabel('Force(pN)')
                        ax.set_xlabel('Time(s)')

                else:
                    RGB_unaltered = stackRGB
                    dx = scanPointer.json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    dt = (scanPointer.timestamps[0,1]-scanPointer.timestamps[0,0]) / 1000000000 #scan time in s
                    totalScanWidth = scanPointer.json['scan volume']['scan axes'][0]['scan width (um)']
                    RGB_altered = modify_rgb_image(RGB_unaltered)
                    #maxTime = len(scanPointer.timestamps[0,:])
                    numberPixels = len(RGB_unaltered)
                    #maxTrueTime = maxTime*dt
                    maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                    if entryScanWidthMax.get() == '-':
                        entryScanWidthMax.delete(0, "end")
                        entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                        labelScanWidthMax = tk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                        labelScanWidthMax.grid(row=5,column=3)
                    if entryYRGBMax.get() == '-':
                        entryYRGBMax.delete(0, "end")
                        entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                        labelYRGBMax = tk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                        labelYRGBMax.grid(row=3,column=3)
                    
                    fig, ax = plt.subplots()
                    
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                    ax.imshow(RGB_altered, **{**default_kwargs})
                    
                    ax.set_ylabel(u'Position(\u03bcm)')
                    ax.set_xlabel(u'Width(\u03bcm)')
                    ax.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                    ax.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])      
                return fig
            
            
            #Function to extract data type that is specific to the Scan data type
            def extractAndPlotScan(h5file,resetBoundsOpt):              
                scanString = typePulldown.get()
                splitScanString = scanString.split('-')
                scanNumber = splitScanString[-1]
                scanPointer = h5file.scans[scanNumber]
                
                reset_axis_items(resetBoundsOpt)
                
                #logical operator to split off if the scan is not one image (aka it is a stack)
                shapeLenRGB = len(scanPointer.rgb_image.shape)
                if shapeLenRGB > 3:
                    fig = extractAndPlotMultipleScans(h5file,scanPointer,resetBoundsOpt)
                    return fig
                
                global labelTimeMax
                global labelDistMax
                global labelYRGBMax
                global labelYForceMax
                global labelScanWidthMax
                
                minTimeScanObj = np.min(scanPointer.timestamps[0,:]).astype(np.int64)
                maxTimeScanObj = np.max(scanPointer.timestamps[0,:]).astype(np.int64)
                
                if plottingOpt.get() == "Both":
                    fig, (ax1, ax2) = plt.subplots(2,constrained_layout=True)
                    RGB_unaltered = scanPointer.rgb_image
                    dx = scanPointer.json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    #dt = (scanPointer.timestamps[0,1]-scanPointer.timestamps[0,0]) / 1000000000 #scan time in s
                    totalScanWidth = RGB_unaltered.shape[1] * dx / 1000
                    RGB_altered = modify_rgb_image(RGB_unaltered)
                    #maxTime = len(scanPointer.timestamps[0,:])
                    numberPixels = len(RGB_unaltered)
                    #maxTrueTime = maxTime*dt
                    maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                    if entryScanWidthMax.get() == '-':
                        entryScanWidthMax.delete(0, "end")
                        entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                        labelScanWidthMax = tk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                        labelScanWidthMax.grid(row=5,column=3)
                    if entryYRGBMax.get() == '-':
                        entryYRGBMax.delete(0, "end")
                        entryYRGBMax.insert(0,str(round(maxTrueDist,3)))
                        labelYRGBMax = tk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                        labelYRGBMax.grid(row=3,column=3)
                    
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                    if grayscaleOpt.get() == "None":
                        ax1.imshow(RGB_altered, **{**default_kwargs})
                    else:
                        ax1.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
                    
                    ax1.set_ylabel(u'Position(\u03bcm)')
                    ax1.set_xlabel(u'Width(\u03bcm)')
                    ax1.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                    ax1.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                    
                    if comboboxForFTvsFD.get() == "Force-Distance":
                        forceData, distData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                        ax2.plot(distData,forceData)
                            
                        if entryDistMax.get() == '-':
                            entryDistMin.delete(0, "end")
                            entryDistMin.insert(0,str(round(min(distData),2)))
                            entryDistMax.delete(0, "end")
                            entryDistMax.insert(0,str(round(max(distData),2)))
                            labelDistMax = tk.Label(frameForAxis,text=str(round(max(distData),2)))
                            labelDistMax.grid(row=2,column=3)
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax2.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                        ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax2.set_ylabel('Force(pN)')
                        ax2.set_xlabel(u'Distance(\u03bcm)')
                        
                        forceString = forceChannelPulldown.get()
                        
                        if len(forceChannelPulldown.get()) > 2:
                            bead_ID = forceString[0]
                            combinedString = 'Bead ' + bead_ID + ' vs. Distance'
                        else:
                            combinedString = 'Channel ' + forceString + ' vs. Distance'
                        ax2.set_title(combinedString)
                        
                    else:
                        forceData, timeData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                        ax2.plot(timeData,forceData)
                        
                        if entryTimeMax.get() == '-':
                            entryTimeMax.delete(0, "end")
                            entryTimeMax.insert(0,str(round(max(timeData),3)))
                            labelTimeMax = tk.Label(frameForAxis,text=str(round(max(timeData),3)))
                            labelTimeMax.grid(row=1,column=3)
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax2.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                        ax2.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax2.set_ylabel('Force(pN)')
                        ax2.set_xlabel('Time(s)')
                        
                        forceString = forceChannelPulldown.get()
                        
                        if len(forceChannelPulldown.get()) > 2:
                            bead_ID = forceString[0]
                            if checkValueDownsampleOpt.get() == 1:
                                combinedString = 'Bead ' + bead_ID + ' Downsampled to ' + entryDownSample.get() + ' Hz'
                            else:
                                combinedString = 'Bead ' + bead_ID
                        else:
                            if checkValueDownsampleOpt.get() == 1:
                                combinedString = forceString + ' Channel Downsampled to ' + entryDownSample.get() + ' Hz'
                            else:
                                combinedString = forceString + ' Channel'
                        ax2.set_title(combinedString)

                        
                elif plottingOpt.get() == "Force Only":
                    fig, ax = plt.subplots(constrained_layout=True)
                    if comboboxForFTvsFD.get() == "Force-Distance":
                        forceData, distData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                        ax.plot(distData,forceData)

                        if entryDistMax.get() == '-':
                            entryDistMin.delete(0, "end")
                            entryDistMin.insert(0,str(round(min(distData),2)))
                            entryDistMax.delete(0, "end")
                            entryDistMax.insert(0,str(round(max(distData),2)))
                            labelDistMax = tk.Label(frameForAxis,text=str(round(max(distData),2)))
                            labelDistMax.grid(row=2,column=3)
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0,"end")
                            entryYForceMin.insert(0,str(round(max(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=3)
                        
                        ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                        ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax.set_ylabel('Force(pN)')
                        ax.set_xlabel(u'Distance(\u03bcm)')
                        
                        forceString = forceChannelPulldown.get()
                        
                    else:
                        forceData, timeData = extract_force_relevant(h5file,timestampsForIndexing=(minTimeScanObj, maxTimeScanObj))
                        ax.plot(timeData,forceData)
                        
                        if entryTimeMax.get() == '-':
                            entryTimeMax.delete(0, "end")
                            entryTimeMax.insert(0,str(round(max(timeData),3)))
                            labelTimeMax = tk.Label(frameForAxis,text=str(round(max(timeData),3)))
                            labelTimeMax.grid(row=1,column=3)
                        
                        if entryYForceMax.get() == '-':
                            entryYForceMin.delete(0, "end")
                            entryYForceMin.insert(0,str(round(min(forceData),2)))
                            entryYForceMax.delete(0, "end")
                            entryYForceMax.insert(0,str(round(max(forceData),2)))
                            labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                            labelYForceMax.grid(row=4,column=4)
                        
                        ax.set_xlim(float(entryTimeMin.get()),float(entryTimeMax.get()))
                        ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                        ax.set_ylabel('Force(pN)')
                        ax.set_xlabel('Time(s)')
                        
                        forceString = forceChannelPulldown.get()

                else:
                    RGB_unaltered = scanPointer.rgb_image
                    dx = scanPointer.json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    #dt = (scanPointer.timestamps[0,1]-scanPointer.timestamps[0,0]) / 1000000000 #scan time in s
                    totalScanWidth = scanPointer.json['scan volume']['scan axes'][0]['scan width (um)']
                    RGB_altered = modify_rgb_image(RGB_unaltered)
                    #maxTime = len(scanPointer.timestamps[0,:])
                    numberPixels = len(RGB_unaltered)
                    #maxTrueTime = maxTime*dt
                    maxTrueDist = dx*numberPixels/1000 #convert to uM
                    
                    
                    if entryScanWidthMax.get() == '-':
                        entryScanWidthMax.delete(0, "end")
                        entryScanWidthMax.insert(0,str(round(totalScanWidth,3)))
                        labelScanWidthMax = tk.Label(frameForAxis,text=str(round(totalScanWidth,3)))
                        labelScanWidthMax.grid(row=5,column=3)
                    
                    if entryYRGBMax.get() == '-':
                        entryYRGBMax.delete(0, "end")
                        entryYRGBMax.insert(0,str(round(maxTrueDist,2)))
                        labelYRGBMax = tk.Label(frameForAxis,text=str(round(maxTrueDist,3)))
                        labelYRGBMax.grid(row=3,column=3)
                    
                    
                    fig, ax = plt.subplots()
                    default_kwargs = dict(extent=[0, totalScanWidth, 0, maxTrueDist], aspect=(RGB_altered.shape[0] / RGB_altered.shape[1]) * (totalScanWidth / maxTrueDist))
                    if grayscaleOpt.get() == "None":
                        ax.imshow(RGB_altered, **{**default_kwargs})
                    else:
                        ax.imshow(RGB_altered, **{**default_kwargs},cmap="gray",norm=matplotlib.colors.NoNorm())
    
                    ax.set_ylabel(u'Position(\u03bcm)')
                    ax.set_xlabel(u'Width(\u03bcm)')
                    ax.set_xlim([float(entryScanWidthMin.get()),float(entryScanWidthMax.get())])
                    ax.set_ylim([float(entryYRGBMin.get()),float(entryYRGBMax.get())])
                return fig

            
            def extractAndPlotFD(h5file,resetBoundsOpt):             
                global labelTimeMax
                global labelDistMax
                global labelYRGBMax
                global labelYForceMax
                global labelScanWidthMax
                
                fdString = typePulldown.get()
                splitFDString = fdString.split('-')
                fdNumber = splitFDString[-1]
                fdPointer = h5file.fdcurves[fdNumber]
                
                reset_axis_items(resetBoundsOpt)
                
                if plottingOpt.get() == "Force Only" or plottingOpt.get() == "Both":
                    if plottingOpt.get() == "Both" or plottingOpt.get() == "RGB Only":
                        print('FD Curves are only supported for Force-Distance plotting only\nIf RGB Data is in the file select the kymo/scan object associated with it to get Force and/or RGB Images')                         
                    
                    fig, ax = plt.subplots(constrained_layout=True)
                    
                    forceString = forceChannelPulldown.get()
                    distanceVarString = whichDistanceValue.get()
                    if len(forceString) < 3:
                        modifiedFDPointer = fdPointer.with_channels(force=forceString,distance=distanceVarString[-1])
                        forceData = modifiedFDPointer.f.data
                        if forceString == "2x":
                            forceData = -1*forceData
                        distData = modifiedFDPointer.d.data
                    else:
                        modifiedFDPointer = fdPointer.with_channels(force=forceString[0],distance=distanceVarString[-1])
                        forceData = modifiedFDPointer.f.data
                        distData = modifiedFDPointer.d.data
                    
                    ax.plot(distData,forceData)
                    
                    if entryDistMax.get() == '-':
                        entryDistMin.delete(0, "end")
                        entryDistMin.insert(0,str(round(min(distData),2)))
                        entryDistMax.delete(0, "end")
                        entryDistMax.insert(0,str(round(max(distData),2)))
                        labelDistMax = tk.Label(frameForAxis,text=str(round(max(distData),2)))
                        labelDistMax.grid(row=2,column=3)
                        
                    if entryYForceMax.get() == '-':
                        entryYForceMin.delete(0, "end")
                        entryYForceMin.insert(0,str(round(min(forceData),2)))
                        entryYForceMax.delete(0, "end")
                        entryYForceMax.insert(0,str(round(max(forceData),2)))
                        labelYForceMax = tk.Label(frameForAxis,text=str(round(max(forceData),2)))
                        labelYForceMax.grid(row=4,column=3)
                        
                    ax.set_xlim(float(entryDistMin.get()),float(entryDistMax.get()))
                    ax.set_ylim(float(entryYForceMin.get()),float(entryYForceMax.get()))
                    ax.set_ylabel('Force(pN)')
                    ax.set_xlabel(u'Distance(\u03bcm)')
                else:
                    print('FD Curves are only supported for Force-Distance plotting only\nIf RGB Data is in the file select the kymo/scan object associated with it to get Force-Time or RGB Images')
                    fig, ax = plt.subplots(constrained_layout=True)
                return fig
            
            #logic to determine if it is the first time this is being plotted
            figureDrawRequest = [directoryPulldown.get(),typePulldown.get(),forceChannelPulldown.get(),whichDistanceValue.get(),comboboxForFTvsFD.get(),checkValueDownsampleOpt.get()]
            
            currentFile = pylake.File(figureDrawRequest[0])
            
            #if resetPlotOpt = 0 then the maximums will not be reset
            resetPlotOpt = 0
            
            #if it is the same file - dont update any exterior parameters
            if figureDrawRequest == lastPlottedFig:
                pass
            #if it is in the same H5 file - dont update metadata just update graph maximums
            elif figureDrawRequest[0] == lastPlottedFig[0]:
                lastPlottedFig[1] = figureDrawRequest[1]
                lastPlottedFig[2] = figureDrawRequest[2]
                lastPlottedFig[3] = figureDrawRequest[3]
                lastPlottedFig[4] = figureDrawRequest[4]
                lastPlottedFig[5] = figureDrawRequest[5]
                resetPlotOpt = 2
            #if it is in a different H5 file - update metadata and update graph maximums
            else:
                lastPlottedFig[0] = figureDrawRequest[0]
                lastPlottedFig[1] = figureDrawRequest[1]
                lastPlottedFig[2] = figureDrawRequest[2]
                lastPlottedFig[3] = figureDrawRequest[3]
                lastPlottedFig[4] = figureDrawRequest[4]
                lastPlottedFig[5] = figureDrawRequest[5]
                resetPlotOpt = 1
                writeMetadata(currentFile)
            
            fileComponent = typePulldown.get()
            splitFileType= fileComponent.split('-')
            filetype = splitFileType[0]
            
            #logical filter to determine what h5 component type is analyzed
            if filetype == "kymos":
                figureReturned = extractAndPlotKymo(currentFile,resetPlotOpt)
            elif filetype == "scans":
                figureReturned = extractAndPlotScan(currentFile,resetPlotOpt)
            else:
                figureReturned = extractAndPlotFD(currentFile,resetPlotOpt)
           
            #Make title of the graph
            figureReturned.suptitle(entryPlotTitle.get(), fontsize=16,va='top')
            #figureReturned.set_size_inches(8,6.4)
            figureReturned.set_dpi(110)
            #make frame, canvas and toolbar objects for the figure
            frameForCanvas = tk.Frame(master,bg="white")
            frameForCanvas.grid(row=0,rowspan=20,column=0,columnspan=2,sticky="nw")
            
            canvas=FigureCanvasTkAgg(figureReturned,master=frameForCanvas)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            canvas.draw()
            
            toolbar = NavigationToolbar2Tk(canvas, frameForCanvas)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            return figureReturned
        
        #Combobox bound function that lists the different .h5 files in that folder, 
        #different file components, and different distance options
        def changeFileComponents(event):
            file = pylake.File(directoryPulldown.get())
            listOfFileTypes = []
            
            for key in file.kymos.keys():
                listOfFileTypes.append('kymos-'+key)
            
            for key in file.scans.keys():
                listOfFileTypes.append('scans-'+key)
            
            for key in file.fdcurves.keys():
                listOfFileTypes.append('fdcurves-'+key)

            typePulldown['values']= listOfFileTypes
            typePulldown.current('0')
            
            #now pull in the distance options and write it to the corresponding combobox pulldown
            distanceGroupPoint = file['Distance']
            listOfDistOptions = []
            for i in distanceGroupPoint:
                listOfDistOptions.append(i)
            whichDistanceValue['values'] = listOfDistOptions
            whichDistanceValue.current('0')
            return
        
        #function to reset file comboboxes upon selecting a new folder        
        def changeH5FileDir(event):
            pointerToDir = filedialog.askdirectory(parent=master, title='Please select a directory to analyze')
        
            os.chdir(pointerToDir)
            h5FileList= glob.glob("*.h5")
            h5FileListTuple = tuple(h5FileList)
            directoryPulldown['values'] = h5FileListTuple
            directoryPulldown.current(0)
            
            #same commands as the changeFileComponents() function
            file = pylake.File(directoryPulldown.get())
            listOfFileTypes = []
            
            for key in file.kymos.keys():
                listOfFileTypes.append('kymos-'+key)
            
            for key in file.scans.keys():
                listOfFileTypes.append('scans-'+key)
            
            for key in file.fdcurves.keys():
                listOfFileTypes.append('fdcurves-'+key)

            typePulldown['values']= listOfFileTypes
            typePulldown.current('0')
            
            #now pull in the distance options
            distanceGroupPoint = file['Distance']
            listOfDistOptions = []
            for i in distanceGroupPoint:
                listOfDistOptions.append(i)
            whichDistanceValue['values'] = listOfDistOptions
            whichDistanceValue.current('0')
        
        #define RGB Constants - change these values based on your normal images
        def define_Global_Defaults():
            global defaultDict
            defaultDict = {}
            defaultDict['red']= 1
            defaultDict['green']= 1
            defaultDict['blue']= 1
            return
        
        #updatePlot button bound event to generate the figure
        def buildPlot(event):
            figureToSave = generateFigure()
            return
        
        #saveImageButton bound event to generate a figure, refresh the image,
        #   and save image/metadata
        def saveFigure(event):
            figureToSave = generateFigure()
            
            imageStringPrefix = (entrySaveFile.get()).replace(" ","_")
            imageSuffix = imageFormatOption.get()
            plt.savefig(imageStringPrefix + '.' +imageSuffix)
            
            #split metadata from the h5 file
            fileName = pylake.File(directoryPulldown.get())
            
            metaData = fileName.description
            
            #open,write, and save metadata File
            metaDataFileString = imageStringPrefix.replace(' ','_')+ '_desc' +'.txt'
            metaDataFile = open(metaDataFileString,'w')
            metaDataFile.write(f'Metadata for image analysis of {directoryPulldown.get()} --- {typePulldown.get()}\n')
            
            metaDataFile.write(metaData)
            
            metaDataFile.write('------------------------------------------\n')
            metaDataFile.write(f'Image Scaling Factors\nR\t{entryRed.get()}\nG\t{entryGreen.get()}\nB\t{entryBlue.get()}\n')
            metaDataFile.write(f'Brightness Shift: {entryBrightness.get()}\n')
            
            arrayForKymoOrScan = (typePulldown.get()).split('-')
            
            if  arrayForKymoOrScan[0] == "scans":
                scanNumber = arrayForKymoOrScan[-1]
                scanPointer = fileName.scans[scanNumber]
                
                if scanPointer.num_frames == 1:
                    timestampArray = scanPointer.timestamps
                    
                    metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0]/1e9} seconds\n")
                    metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[:,:]) - timestampArray[0,0])/1e9} seconds\n")
                else:
                    timestampArray = scanPointer.timestamps
                    metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0,0]/1e12} seconds\n")
                    metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[0,:,:]) - timestampArray[0,0,0])/1e9} seconds\n")

            return
        
        """
        This function has been borrowed from Michael Wasserman's code to extract the image data
        for the current file being analyzed. This code is different from the batch extract code (h5_image_extract.py)
        is that no force vs. RGB summary plot are made. This functionality is still availiable using the save plot
        button to get force-correlated RGB data (or force only data)
        This makes the image extract button much faster because you do not have to read in large force datasets
        """
        def extractImageCTrap(event):
            temp_file_name = directoryPulldown.get()
            
            # extract and save experimental description
            def save_exp_desc(exp_desc,filename_without_extension):
                imageStringPrefix = filename_without_extension
                
                fileName = pylake.File(directoryPulldown.get())
                
                metaDataFileString = imageStringPrefix.replace(' ','_')+ '_desc' +'.txt'
                metaDataFile = open(metaDataFileString,'w')
                metaDataFile.write(f'Metadata for image analysis of {directoryPulldown.get()} --- {typePulldown.get()}\n')
                
                metaDataFile.write(exp_desc)
            
                metaDataFile.write('------------------------------------------\n')
                metaDataFile.write(f'Image Scaling Factors\nR\t{entryRed.get()}\nG\t{entryGreen.get()}\nB\t{entryBlue.get()}\n')
                metaDataFile.write(f'Brightness Shift: {entryBrightness.get()}\n')
            
                arrayForKymoOrScan = (typePulldown.get()).split('-')
            
                if  arrayForKymoOrScan[0] == "scans":
                    scanNumber = arrayForKymoOrScan[-1]
                    scanPointer = fileName.scans[scanNumber]
                
                    if scanPointer.num_frames == 1:
                        timestampArray = scanPointer.timestamps
                        metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0]/1e9} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[:,:]) - timestampArray[0,0])/1e9} seconds\n")
                    else:
                        timestampArray = scanPointer.timestamps
                        metaDataFile.write(f"First Timestamp Value: {timestampArray[0,0,0]/1e9} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {(np.max(timestampArray[0,:,:]) - timestampArray[0,0,0])/1e9} seconds\n")
            
            # save scans without labeled axes. Important for conserving pixel num etc
            def save_image_for_CTrapViewer(image_num, image_type):
                if image_type == 'kymo':
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.kymos[s].json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    #dt = tempfile.kymos[s].json['scan volume']['scan axes'][0]['scan time (ms)'] #scan time in ms
                    # Note: the above dt extraction is incorrect if you manually add a wait time between line scans
                    dt = (tempfile.kymos[s].timestamps[0,1]-tempfile.kymos[s].timestamps[0,0])/1e6 #scan time in ms
                    # save kymograph without labeled axes. Important for conserving pixel num etc
                    filename_png = filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms" + ".png"
                    tempfile.kymos[s].save_tiff(filename_png)
                elif image_type == 'scan':
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.scans[t].json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    if tempfile.scans[t].num_frames == 1: #single 2D scan
                        dt = (tempfile.scans[t].timestamps[0,1]-tempfile.scans[t].timestamps[0,0])/1e6 #line time in ms
                        image_slice = tempfile.scans[t].rgb_image
                        image_slice_mod = modify_rgb_image(image_slice)
                        heightScan = image_slice.shape[1] * dx / 1e3
                        totalScanWidth =  tempfile.scans[t].json['scan volume']['scan axes'][0]['scan width (um)']
                        default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_mod.shape[0] / image_slice_mod.shape[1]) * (totalScanWidth/ heightScan))
                            
                        filename_png = filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms" + ".png"
                        fig, ax = plt.subplots(constrained_layout=True)
                        ax.imshow(image_slice_mod, **{**default_kwargs})
                        plt.axis('off')
                        plt.savefig(filename_png,bbox_inches="tight", pad_inches = 0)
                        
                    elif tempfile.scans[t].num_frames > 1: #image stack
                        dt = (tempfile.scans[t].timestamps[0,1,0]-tempfile.scans[t].timestamps[0,0,0])/1e6 #line time in ms

                        for i in range(tempfile.scans[t].num_frames):
                            frame_num = i + 1
                            filename_png = filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms f" + str(frame_num) + ".png"
                            image_slice = tempfile.scans[t].rgb_image[i,:,:,:]
                            heightScan = image_slice.shape[1] * dx / 1e3
                            image_slice_int = modify_rgb_image(image_slice) #convert to uint8 for saving as .png <- JWW changed this, might need to change back
                            totalScanWidth =  tempfile.scans[t].json['scan volume']['scan axes'][0]['scan width (um)']
                            
                            #listOfInterpols = ['none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric','catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
                            default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_int.shape[0] / image_slice_int.shape[1]) * (totalScanWidth/ heightScan),interpolation="nearest")
                            
                            fig, ax = plt.subplots(constrained_layout=True)
                            ax.imshow(image_slice_int, **{**default_kwargs})
                            plt.axis('off')
                            #filename_png = item + filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms f" + str(frame_num) + ".png"
                            plt.savefig(filename_png,bbox_inches = 'tight', pad_inches = 0)
          
            tempfile = pylake.File(temp_file_name) # load h5 file
            filename_without_extension = tempfile.h5.filename.replace(".h5", "")  # "file.h5" -> "file"
            exp_desc = tempfile.description        
            save_exp_desc(exp_desc, filename_without_extension) # save .txt with experimental description
                
            # loop to extract each kymo 
            for s in list(tempfile.kymos): 
                image_type = 'kymo'                                 
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_CTrapViewer(s, image_type)
               
            # loop to extract each 2D scan 
            for t in list(tempfile.scans): 
                image_type = 'scan'            
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_CTrapViewer(t, image_type)
            return
        
        def extractImageImageJ(event):
            temp_file_name = directoryPulldown.get()
            
            # extract and save experimental description
            def save_exp_desc(exp_desc,filename_without_extension):
                imageStringPrefix = filename_without_extension
                
                fileName = pylake.File(directoryPulldown.get())
                
                metaDataFileString = imageStringPrefix.replace(' ','_')+ '_desc' +'.txt'
                metaDataFile = open(metaDataFileString,'w')
                metaDataFile.write(f'Metadata for image analysis of {directoryPulldown.get()} --- {typePulldown.get()}\n')
                
                metaDataFile.write(exp_desc)
            
                metaDataFile.write('------------------------------------------\n')
                metaDataFile.write(f'Image Scaling Factors\nR\t{entryRed.get()}\nG\t{entryGreen.get()}\nB\t{entryBlue.get()}\n')
                metaDataFile.write(f'Brightness Shift: {entryBrightness.get()}\n')
            
                arrayForKymoOrScan = (typePulldown.get()).split('-')
                if  arrayForKymoOrScan[0] == "scans":
                    scanNumber = arrayForKymoOrScan[-1]
                    scanPointer = fileName.scans[scanNumber]
                
                    
                    if scanPointer.num_frames == 1:
                        timestampArray = scanPointer.timestamps
                        firstTimestamp = timestampArray[0,0]/1e9
                        timestepBetweenImages = (np.max(timestampArray[:,:]) - timestampArray[0,0])/1e9
                        metaDataFile.write(f"First Timestamp Value: {firstTimestamp} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {timestepBetweenImages} seconds\n")
                    else:
                        timestampArray = scanPointer.timestamps
                        firstTimestamp = timestampArray[0,0,0]/1e9
                        timestepBetweenImages = (np.max(timestampArray[0,:,:]) - timestampArray[0,0,0])/1e9
                        metaDataFile.write(f"First Timestamp Value: {firstTimestamp} seconds\n")
                        metaDataFile.write(f"Image Acquisition Time: {timestepBetweenImages} seconds\n")
                
                return
            
            # save scans without labeled axes. Important for conserving pixel num etc
            def save_image_for_ImageJ(image_num, image_type):
                if image_type == 'kymo':
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.kymos[s].json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm
                    #dt = tempfile.kymos[s].json['scan volume']['scan axes'][0]['scan time (ms)'] #scan time in ms
                    # Note: the above dt extraction is incorrect if you manually add a wait time between line scans
                    dt = (tempfile.kymos[s].timestamps[0,1]-tempfile.kymos[s].timestamps[0,0])/1e6 #scan time in ms
                    # save kymograph without labeled axes. Important for conserving pixel num etc
                    filename_png = filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms" + ".png"
                    tempfile.kymos[s].save_tiff(filename_png)
                elif image_type == 'scan':
                    # Add dx and dt to filename for CTrapViewer purposes
                    dx = tempfile.scans[t].json['scan volume']['scan axes'][0]['pixel size (nm)'] #pixel size in nm

                    if tempfile.scans[t].num_frames == 1: #single 2D scan
                        dt = (tempfile.scans[t].timestamps[0,1]-tempfile.scans[t].timestamps[0,0])/1e6 #line time in ms
                        image_slice = tempfile.scans[t].rgb_image
                        image_slice_mod = modify_rgb_image(image_slice)
                        heightScan = image_slice.shape[0] * dx / 1e3
                        #totalScanWidth =  tempfile.scans[t].json['scan volume']['scan axes'][0]['scan width (um)']
                        totalScanWidth = image_slice.shape[1] * dx / 1e3
                        default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_mod.shape[0] / image_slice_mod.shape[1]) * (totalScanWidth/ heightScan))
                            
                        filename_png = filename_without_extension + " ImageJ " +" dx " + str(dx) + "nm dt " + str(dt) + "ms" + ".png"
                        fig, ax = plt.subplots(constrained_layout=True)
                        ax.imshow(image_slice_mod, **{**default_kwargs})
                        ax.set_ylabel(u'Position(\u03bcm)')
                        ax.set_xlabel(u'Scan Width(\u03bcm)')
                        fig.suptitle(f' ', y = -0.04)
                        plt.savefig(filename_png,bbox_inches="tight", pad_inches = 0)
                        
                    elif tempfile.scans[t].num_frames > 1: #image stack
                        dt = (tempfile.scans[t].timestamps[0,1,0]-tempfile.scans[t].timestamps[0,0,0])/1e6 #line time in ms
                        
                        initialTimeTotal = tempfile.scans[t].timestamps[0,0,0]
                        
                        for i in range(tempfile.scans[t].num_frames):
                            frame_num = i + 1
                            filename_png = filename_without_extension + " dx " + str(dx) + "nm dt " + str(dt) + "ms f" + str(frame_num) + ".png"
                            image_slice = tempfile.scans[t].rgb_image[i,:,:,:]
                            initialTimeScan = tempfile.scans[t].timestamps[i,0,0]
                            heightScan = image_slice.shape[0] * dx / 1e3
                            image_slice_int = modify_rgb_image(image_slice) #convert to uint8 for saving as .png <- JWW changed this, might need to change back
                            #totalScanWidth =  tempfile.scans[t].json['scan volume']['scan axes'][0]['scan width (um)']
                            totalScanWidth = image_slice.shape[1] * dx / 1e3
                            #listOfInterpols = ['none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric','catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
                            default_kwargs = dict(extent=[0, totalScanWidth, 0, heightScan], aspect=(image_slice_int.shape[0] / image_slice_int.shape[1]) * (totalScanWidth/ heightScan),interpolation="nearest")
                            
                            fig, ax = plt.subplots(constrained_layout=True)
                            ax.imshow(image_slice_int, **{**default_kwargs})
                            ax.set_ylabel(u'Position(\u03bcm)')
                            ax.set_xlabel(u'Scan Width(\u03bcm)')
                            fig.suptitle(f'{round((initialTimeScan - initialTimeTotal) / 1e9,2)} s', y = -0.04)
                            plt.savefig(filename_png,bbox_inches = 'tight', pad_inches = 0)
          
            tempfile = pylake.File(temp_file_name) # load h5 file
            filename_without_extension = tempfile.h5.filename.replace(".h5", "")  # "file.h5" -> "file"
            exp_desc = tempfile.description        
            save_exp_desc(exp_desc, filename_without_extension) # save .txt with experimental description
                
            # loop to extract each kymo 
            for s in list(tempfile.kymos): 
                image_type = 'kymo'                                 
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_ImageJ(s, image_type)
               
            # loop to extract each 2D scan 
            for t in list(tempfile.scans): 
                image_type = 'scan'            
                # save kymograph as png for viewing in CTrapViewer
                save_image_for_ImageJ(t, image_type)
            return
        
        """
        This function has been borrowed from Michael Wasserman's code to extract the desired force data
        for the current file being analyzed.
        This task is computationally very expensive if you want the HF data - be wary of this
        """
        def extractForceForOrigin(event):
            def extractForceCommand(settings_list):
                # Extract force data
                filename = directoryPulldown.get()
                temp_file = pylake.File(filename)
                filename_no_extension = filename.replace(".h5",'')
                force1xHF = temp_file['Force HF']['Force 1x']
                force1yHF = temp_file['Force HF']['Force 1y']
                force2xHF = temp_file['Force HF']['Force 2x']
                force2yHF = temp_file['Force HF']['Force 2y']
                sample_rate = force1xHF.sample_rate # Hz
                downsampled_rate = float(entryDownSample.get())
                
                force1x_downsamp = force1xHF.downsampled_by(int(sample_rate/downsampled_rate))
                force1y_downsamp = force1yHF.downsampled_by(int(sample_rate/downsampled_rate))
                force2x_downsamp = force2xHF.downsampled_by(int(sample_rate/downsampled_rate))
                force2y_downsamp = force2yHF.downsampled_by(int(sample_rate/downsampled_rate))
                pooled_force_data = [force1xHF, force1yHF, force2xHF, force2yHF, force1x_downsamp,
                                     force1y_downsamp, force2x_downsamp, force2y_downsamp]
                
                #extract HF and downsampled time stamps
                time = force1xHF.timestamps/1e9 # time traces (seconds)
                time = time - time[0]
                time_downsamp = force1x_downsamp.timestamps/1e9
                time_downsamp = time_downsamp - time_downsamp[0]

                #indices for what force channels to take                
                force_ind = [i for i, v in enumerate(settings_list) if v[:][1] == True]
                if len(force_ind) > 0:
                    all_column_names, force_bool = zip(*settings_list) #split pairs
                    column_names =  []
                    forces_to_save = []
                if any(force_bool[0:3]): # if any HF data is to be saved
                    forces_to_save.append(time)
                    column_names.append('Time (s)')
                    for i in range(4):
                        if force_bool[i] == True:
                            forces_to_save.append(pooled_force_data[i].data)
                            column_names.append(all_column_names[i])
                if any(force_bool[4:7]): # if any downsampled data is to be saved
                    forces_to_save.append(time_downsamp)
                    column_names.append('Time downsamp (s)')               
                    for j in range(4,8):
                        if force_bool[j] == True:
                            forces_to_save.append(pooled_force_data[j].data)
                            column_names.append(all_column_names[j]) 
                else:
                    print('No data channels were selected')
                    return
                
                #save metadata to a .txt file
                exp_desc = temp_file.description
                filename_meta = filename_no_extension + "_desc.txt"
                if len(exp_desc) > 0 and not os.path.exists(filename_meta):
                    meta_file = open(filename_meta, "w")
                    meta_file.write(exp_desc)
                    meta_file.close()
                
                ## preallocate pd DataFrame with #indices = max length of force data
                grouped_data = pd.DataFrame(index=range(len(forces_to_save[0])))
                for i in range(len(forces_to_save)):
                    # if shorter than index length, append NaN: 
                    if len(forces_to_save[i]) < len(forces_to_save[0]):
                        forces_to_save[i] = np.append(forces_to_save[i], np.repeat(np.nan, len(forces_to_save[0])-len(forces_to_save[i])))
                    #transpose force arrays for saving in columns
                    forces_to_save[i]=forces_to_save[i].reshape(-1,1)
                    grouped_data.insert(i, column_names[i], forces_to_save[i])
                
                grouped_data = grouped_data.fillna('') # fill NaN with ''
                filename_force = filename_no_extension + "_force.csv"
                grouped_data.to_csv(filename_force, sep=',', index=False, header=True)
                return
            
            def extractAndDestroyForceOrigin():
                settings_list = [('Force1x HF', force1xHF_var.get()), 
                     ('Force1y HF', force1yHF_var.get()), 
                     ('Force2x HF', force2xHF_var.get()), 
                     ('Force2y HF', force2yHF_var.get()),                      
                     ('Force1x downsamp', force1x_downsamp_var.get()), 
                     ('Force1y downsamp', force1y_downsamp_var.get()), 
                     ('Force2x downsamp', force2x_downsamp_var.get()), 
                     ('Force2y downsamp',force2y_downsamp_var.get())]
                extractForceCommand(settings_list)
                forceRoot.destroy()
                return
            
            #Pop Up Menu
            forceRoot = tk.Toplevel()
            forceRoot.wm_title("Force Extract Option")
            
            # Initialize variables for save option checkboxes
            force1xHF_var = tk.BooleanVar()
            force1yHF_var = tk.BooleanVar()
            force2xHF_var = tk.BooleanVar()
            force2yHF_var = tk.BooleanVar()    
            force1x_downsamp_var = tk.BooleanVar()
            force1x_downsamp_var.set(True) # most common save option, so selected by default
            force1y_downsamp_var = tk.BooleanVar()
            force2x_downsamp_var = tk.BooleanVar()
            force2y_downsamp_var = tk.BooleanVar()
            
            # Establish checkboxes for each save option
            tk.Checkbutton(forceRoot, text="Force1x_HF",variable=force1xHF_var).grid(row=1, sticky="w")
            tk.Checkbutton(forceRoot, text="Force1y_HF",variable=force1yHF_var).grid(row=2, sticky="w")
            tk.Checkbutton(forceRoot, text="Force2x_HF",variable=force2xHF_var).grid(row=3, sticky="w")
            tk.Checkbutton(forceRoot, text="Force2y_HF",variable=force2yHF_var).grid(row=4, sticky="w")
            tk.Checkbutton(forceRoot, text="Force1x_downsamp",variable=force1x_downsamp_var).grid(row=5, sticky="w")
            tk.Checkbutton(forceRoot, text="Force1y_downsamp",variable=force1y_downsamp_var).grid(row=6, sticky="w")
            tk.Checkbutton(forceRoot, text="Force2x_downsamp",variable=force2x_downsamp_var).grid(row=7, sticky="w")
            tk.Checkbutton(forceRoot, text="Force2y_downsamp",variable=force2y_downsamp_var).grid(row=8, sticky="w")
            tk.Button(forceRoot, text='Extract data', command=extractAndDestroyForceOrigin).grid(row=10, sticky="w", pady=8)
            return
        
        #showImageHistogram bound button event to display the RGB histogram
        def showImageHistogram(event):
            #destroyHistImage bound event to quit the popup menu
            def destroyHistImage(event):
                histRoot.destroy()
                return
            
            inputComponent = typePulldown.get()
            kymoOrScan = inputComponent.split('-')          
            h5file = pylake.File(directoryPulldown.get())
            histRoot = tk.Toplevel(bg="white")
            histRoot.wm_title("RGB Photon Count Distribution Post Modifications")
            
            #configure button to exit the popup menu
            quitHistButton = tk.Button(histRoot, text="Quit\nHistogram\nViewer",width=10,height=5)
            quitHistButton.grid(row=4,column=1)
            quitHistButton.bind("<Button-1>",destroyHistImage)
            
            if kymoOrScan[0] == 'kymos':
                rgb_data = h5file.kymos[kymoOrScan[-1]].rgb_image
            elif kymoOrScan[0] =='scans':
                rgb_data = h5file.scans[kymoOrScan[-1]].rgb_image
            else:
                print('No RGB data found in FD files!')
            
            #run the RGB image modifier
            rgb_image_modified = modify_rgb_image(rgb_data)
            
            #configure histogram plot
            redData = rgb_image_modified[:,:,0].flatten()
            greenData = rgb_image_modified[:,:,1].flatten()
            blueData = rgb_image_modified[:,:,2].flatten()
            
            channel_max_vals = np.array([np.max(redData),np.max(greenData),np.max(blueData)])
            
            figHist, ax = plt.subplots()
            
            maxVal = np.max(channel_max_vals)
            
            numBins = int(maxVal/5)
            
            binRange = (0,maxVal)
            
            ax.hist(redData, bins=numBins, range = binRange, density=False,color="red",label="Red",alpha=0.5, histtype='step')
            ax.hist(greenData, bins=numBins, range = binRange, density=False,color="green",label="Green",alpha=0.5, histtype='step')
            ax.hist(blueData, bins=numBins, range = binRange, density=False,color="blue",label="Blue",alpha=0.5, histtype='step')
            
            ax.legend(loc="upper right") 
            
            figHist.suptitle("RGB Photon Counts")
            ax.set_ylabel('Occurences')
            ax.set_xlabel('Photon Counts')
            
            #write the canvas figure
            frameForCanvas = tk.Frame(histRoot,bg="white")
            frameForCanvas.grid(row=0,column=0,rowspan=9)
            canvas=FigureCanvasTkAgg(figHist,master=frameForCanvas)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            canvas.draw()
            return
        
        #quitButton bound event to destory the tkinter window and exit out of python
        def totalQuit(event):
            master.destroy()
            quit()
            return
        
        self.master = master
        #build the simple GUI
        define_Global_Defaults()
        ###NOTE: A lot of these entries/comboboxes/ are defined as globals 
        #add directory system - values to be assigned dynamically later
        frameForFileAccess = tk.Frame(master)
        frameForFileAccess.grid(row=0,column=2,columnspan=2,rowspan=2,sticky='w',padx=10)
        tk.Label(frameForFileAccess,text="H5 Files in Directory",height=1,font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2,sticky='w')
        
        directoryPulldown = tk.ttk.Combobox(frameForFileAccess,width=48)
        directoryPulldown.grid(row=1,column=0,columnspan=4)
        tk.Label(frameForFileAccess,text="File Components",font=('Helvetica', 10, 'bold')).grid(row=2,column=0,columnspan=1,sticky='w')
        
        typePulldown = tk.ttk.Combobox(frameForFileAccess,width=15)
        typePulldown.grid(row=3,column=0,columnspan=1,sticky='w')
            
        ###add color scaling options
        #Labels
        frameForColorOpt = tk.Frame(master)
        frameForColorOpt.grid(row=4,column=3,rowspan=3,columnspan=1,sticky='w',padx=10)
        tk.Label(frameForColorOpt,text="Photon Count Multiplier",font=('Helvetica', 10, 'bold')).grid(row=0 ,column=0,columnspan=2,sticky='nw')
        #RGBscalingFrame = tk.Frame(frameForColorOpt)
        #RGBscalingFrame.grid(row=1,column=1)
        tk.Label(frameForColorOpt,text="Red Scalar").grid(row=1, column=0,sticky='w')
        tk.Label(frameForColorOpt,text="Green Scalar").grid(row=2, column=0,sticky='w')
        tk.Label(frameForColorOpt,text="Blue Scalar").grid(row=3, column=0,sticky='w')
        tk.Label(frameForColorOpt,text="Brightness").grid(row=4, column=0,sticky='w')
            
        #Entries for color analysis
        entryRed = tk.Entry(frameForColorOpt,width=5)
        entryRed.insert(0,str(defaultDict['red']))
        entryRed.grid(row=1,column=1)
        entryGreen = tk.Entry(frameForColorOpt,width=5)
        entryGreen.insert(0,str(defaultDict['green']))
        entryGreen.grid(row=2,column=1)
        entryBlue = tk.Entry(frameForColorOpt,width=5)
        entryBlue.insert(0,str(defaultDict['blue']))
        entryBlue.grid(row=3,column=1)
        entryBrightness = tk.Entry(frameForColorOpt,width=5)
        entryBrightness.grid(row=4,column=1)
        entryBrightness.insert(0,'0')
        tk.Label(frameForColorOpt,text="Greyscale Options").grid(row=5 ,column=0,columnspan=1,sticky='nw')
            
        grayscaleOpt = tk.ttk.Combobox(frameForColorOpt,values=['None','Red','Green','Blue'],width=6)
        grayscaleOpt.current('0')
        grayscaleOpt.grid(row=5, column=1,sticky='w')
            
        bugFixButton= tk.BooleanVar()
        tk.Checkbutton(frameForColorOpt, text="Fix Image Reconstruction?\n(Pylake <v0.6.0 Bug)",variable=bugFixButton).grid(row=6, column=0, columnspan=2)
            
        #placeholder for future image analysis tool
            
        ###add force options
        frameForForceOptions = tk.Frame(master)
        frameForForceOptions.grid(row=2,column=3,rowspan=2,columnspan=1,sticky='w',pady=10,padx=10)
        tk.Label(frameForForceOptions,text="Force Options",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2,sticky='w')
        tk.Label(frameForForceOptions,text="Sampling(Hz)",).grid(row=1, column=0,sticky='w')
        tk.Label(frameForForceOptions,text="Channel",).grid(row=2, column=0,sticky='w')
        
        entryDownSample = tk.Entry(frameForForceOptions,width=4)
        entryDownSample.grid(row=1,column=1)
        entryDownSample.insert(0,'100')
        forceChannelPulldown = tk.ttk.Combobox(frameForForceOptions,values=['1x','1y','2x','2y','1-Both','2-Both'],width=4)
        forceChannelPulldown.grid(row=2,column=1)
        forceChannelPulldown.current('2')
        checkValueDownsampleOpt = tk.IntVar(value=1)
        forceOptions = tk.Checkbutton(frameForForceOptions,text="Downsample?",var=checkValueDownsampleOpt)
        forceOptions.grid(row=3,column=0,columnspan=2)
            
        ###add plotting options
        frameForPlotting = tk.Frame(master)
        frameForPlotting.grid(row=2,column=2,rowspan=2,columnspan=1,sticky='w',pady=10)
        tk.Label(frameForPlotting,text="Plotting Options",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2,sticky='w')
        tk.Label(frameForPlotting,text="Which Plot(s)").grid(row=1,column=0,sticky='w')
        tk.Label(frameForPlotting,text="Force Plot").grid(row=2,column=0,sticky='w')
        tk.Label(frameForPlotting,text="Distance Data").grid(row=3,column=0,sticky='w')
        plottingOpt = tk.ttk.Combobox(frameForPlotting,values=['Both','RGB Only','Force Only'],width=10)
        plottingOpt.grid(row=1,column=2)
        plottingOpt.current('1')
        comboboxForFTvsFD = tk.ttk.Combobox(frameForPlotting,values=['Force-Time','Force-Distance'],width=10)
        comboboxForFTvsFD.grid(row=2,column=2)
        comboboxForFTvsFD.current('0')
        whichDistanceValue = tk.ttk.Combobox(frameForPlotting,width=10)
        whichDistanceValue.grid(row=3,column=2)
        
        ###axis scaling interface
        frameForAxis = tk.Frame(master)
        frameForAxis.grid(row=4,column=2,rowspan=3,columnspan=1,sticky='w')
        tk.Label(frameForAxis,text="Axis Options",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=4,sticky='w')
        tk.Label(frameForAxis,text="Time").grid(row=1,column=0,sticky='w')
        tk.Label(frameForAxis,text="Distance").grid(row=2,column=0,sticky='w')
        tk.Label(frameForAxis,text="Y RGB").grid(row=3,column=0,sticky='w')
        tk.Label(frameForAxis,text="Y Force").grid(row=4,column=0,sticky='w')
        tk.Label(frameForAxis,text="Scan Width").grid(row=5,column=0,sticky='w')
            
        axisEntryWidth = 6
        padAxisEntry = 2
        #Minimum axis entries - set by default to zero
        entryTimeMin = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryTimeMin.grid(row=1,column=1,padx=padAxisEntry)
        entryTimeMin.insert(0,'0')
        entryDistMin = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryDistMin.grid(row=2,column=1,padx=padAxisEntry)
        entryDistMin.insert(0,'0')
        entryYRGBMin = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryYRGBMin.grid(row=3,column=1,padx=padAxisEntry)
        entryYRGBMin.insert(0,'0')
        entryYForceMin = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryYForceMin.grid(row=4, column=1,padx=padAxisEntry)
        entryYForceMin.insert(0,'0')
        entryScanWidthMin = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryScanWidthMin.grid(row=5,column=1,padx=padAxisEntry)
        entryScanWidthMin.insert(0,'0')
            
        #Maximum axis entries - not defined by any default since it depends on the file
        entryTimeMax = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryTimeMax.grid(row=1,column=2,padx=padAxisEntry)
        entryTimeMax.insert(0,'-')
        entryDistMax = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryDistMax.grid(row=2,column=2,padx=padAxisEntry)
        entryDistMax.insert(0,'-')
        entryYRGBMax = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryYRGBMax.grid(row=3,column=2,padx=padAxisEntry)
        entryYRGBMax.insert(0,'-')
        entryYForceMax = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryYForceMax.grid(row=4, column=2,padx=padAxisEntry)
        entryYForceMax.insert(0,'-')
        entryScanWidthMax = tk.Entry(frameForAxis,width=axisEntryWidth)
        entryScanWidthMax.grid(row=5,column=2,padx=padAxisEntry)
        entryScanWidthMax.insert(0,'-')
        
        #Labels and entries for Plot Title and Image/metadata
        tk.Label(frameForFileAccess,text="Plot Title",font=('Helvetica', 10, 'bold')).grid(row=2,column=1,sticky='w')
        tk.Label(frameForFileAccess,text="File Name",font=('Helvetica', 10, 'bold')).grid(row=3,column=1,sticky='w')
        entryPlotTitle = tk.Entry(frameForFileAccess)
        entryPlotTitle.grid(row=2,column=2)
        #entryPlotTitle.insert(0,)    
        entrySaveFile = tk.Entry(frameForFileAccess)
        entrySaveFile.grid(row=3,column=2)
            
        frameForImageFormat = tk.Frame(frameForFileAccess)
        frameForImageFormat.grid(row=4,column=1,columnspan=2)
        tk.Label(frameForImageFormat,text="Image Format").pack(side="left")
        imageFormatOption = tk.ttk.Combobox(frameForImageFormat,values=['PNG','PDF','TIFF','JPEG','SVG'],width=5)
        imageFormatOption.pack(side="right")
        imageFormatOption.current('0')
            
        buttonHeight = 3
        buttonWidth = 15
        #Buttons to link to commands
        buttonFrame = tk.Frame(master)
        buttonFrame.grid(row=0,column=4,rowspan=9)
        updateH5File= tk.Button(frameForFileAccess,text="Change Folder")
        updateH5File.grid(row=0,column=2,columnspan=2,pady=1)
        updatePlot = tk.Button(buttonFrame,text="Draw Plot",width=buttonWidth,height=buttonHeight)
        updatePlot.pack(side="top",padx=4,pady=2)
        exportForceButton = tk.Button(buttonFrame,text="Export Force\nFor Origin",width=buttonWidth,height=buttonHeight)
        exportForceButton.pack(side="top",padx=4,pady=2)
        exportImageButton = tk.Button(buttonFrame,text="Export Image For\n CTrapViewer",width=buttonWidth,height=buttonHeight)
        exportImageButton.pack(side="top",padx=4,pady=2)
        exportForImageJ = tk.Button(buttonFrame,text="Export Image For\n ImageJ",width=buttonWidth,height=buttonHeight) #button to quit Tkinter GUI
        exportForImageJ.pack(side="top",padx=4,pady=2)
        showRGBHist = tk.Button(buttonFrame,text="Show RGB\nHistogram",width=buttonWidth,height=buttonHeight)
        showRGBHist.pack(side="top",padx=4,pady=2)
        saveImageButton = tk.Button(buttonFrame,text="Save GUI Image",width=buttonWidth,height=buttonHeight)
        saveImageButton.pack(side="top",padx=4,pady=2)
        quitButton = tk.Button(buttonFrame,text="Quit?",width=buttonWidth,height=buttonHeight) #button to quit Tkinter GUI
        quitButton.pack(side="top",padx=4,pady=2)
        
        #Inputs and Labels for metadata
        frameForMetadata = tk.Frame(master)
        frameForMetadata.grid(row=7,column=2,rowspan=3,columnspan=2,sticky='w')
        tk.Label(frameForMetadata,text="Metadata",font=('Helvetica', 10, 'bold')).grid(row=0,column=0,columnspan=2,sticky='w')
        
        lastPlottedFig = ['','','','','','']
            
        frameForSlider = tk.Frame(master)
        frameForSlider.grid(row=10,column=2,columnspan=1,sticky="w")
        tk.Label(frameForSlider,text="Scan Image Frame",font=('Helvetica', 10, 'bold'),justify="left").grid(row=0,column=0,sticky="w")
        scaleForStack = tk.Scale(frameForSlider,from_=0,to=100,orient=tk.HORIZONTAL)
        scaleForStack.grid(row=1,column=0)
        multiScanPlotOpt = tk.IntVar(value=1)
        plotScanTime = tk.Checkbutton(frameForSlider,text="Highlight Scan Range?",var=multiScanPlotOpt,anchor="w",justify=tk.LEFT)
        plotScanTime.select()
        plotScanTime.grid(row=2,column=0,sticky="w")
                
        #bind the proper buttons to the desired functions
        quitButton.bind("<ButtonRelease-1>",totalQuit)
        updateH5File.bind("<ButtonRelease-1>",changeH5FileDir)
        exportForceButton.bind("<ButtonRelease-1>",extractForceForOrigin)
        exportImageButton.bind("<ButtonRelease-1>",extractImageCTrap)
        exportForImageJ.bind("<ButtonRelease-1>",extractImageImageJ)
        directoryPulldown.bind("<<ComboboxSelected>>",changeFileComponents)
        updatePlot.bind("<ButtonRelease-1>",buildPlot)
        master.bind("<Return>",buildPlot)
        showRGBHist.bind("<ButtonRelease-1>", showImageHistogram)
        saveImageButton.bind("<ButtonRelease-1>",saveFigure)
        master.bind("<Escape>",totalQuit)
    

#call commands to open the gui class and loop continuosly
root = tk.Tk()
root.title('C-TrapVis v1.0.0')
my_gui = CTrapGUI(root) #call the Application class
root.mainloop()