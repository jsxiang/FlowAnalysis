# FlowAnalysis

Flow cytomery analysis scripts used to analyze exported .fcs files. Data are generated on the MACSQuant VYB analyzer. fca_fcsread.m is required for extracting data from .fcs files (see fca_fcsread_license.txt). 

Input files indicating experiment conditions are to be prepared as a text file. See input1mM.txt and input5mM.txt as example. 

analysis.m is the main analysis file that uses analyzeflow.m to find relative fluorescence values as mCherry/BFP, or GFP/mCherry of every sample, and then uses getVYBdata.m to organize the fluorescence values by determining the mean and standard deviation for each sequence tested. 
