#filter sensor data

%reset -f 

import os
import pandas as pd
import numpy as np

#Create butterworth lowpass filter. From Scipy
from scipy import signal
fs=100       #sampling rate in hz
fc=20        #cutoff frequency of filter in hz
order=4
Wn=fc/(fs/2) 
#critical frequencies for butterworth filter (the frequency where the magnitude response of the filter is 1/sqrt(2)), expressed as the fraction of the Nyquist frequency, which is half the sampling frequency. For digital filters, the cutoff frequency must lie between 0 and 1, where 1 corresponds to the nyquist rate)
#For digital filters, Wn are in the same units as fs. By default, fs is 2 half-cycles/sample, so these are normalized from 0 to 1, where 1 is the Nyquist frequency. (Wn is thus in half-cycles / sample.)

bf,af=signal.butter(order, Wn, btype='low', analog=False)

#Loop through data

Phases=["P", "SF1", "R1", "R2", "SF2"]

DataDirectory=r"XXXXX"
OutputDirectory=r"XXXXX"

#loop across all datafiles


for Phase in Phases:
    
    os.chdir(DataDirectory+Phase)
    Datafiles=os.listdir(DataDirectory+Phase)

    for f in Datafiles:

        #skip invalid filenames
        if f.startswith("._"): continue
        elif f.endswith(".csv")==False: continue

        #obtain ID information - in first 4 characters of filename
        ID=f[0:4]
            
        #if ID=="XXXXX":

        #obtain suffix information
        #Specify file suffix information
        suffix= ["Wrist-Right",
             "Wrist-Left",
             "Ankle-Right",
             "Ankle-Left",
             "Torso"]

        for suf in suffix:
            if suf in f: break

        #Load data
        data=pd.read_csv(f, engine="python")

        #filter x y z vector
        data["FreeAcc_X_20hz"]=signal.filtfilt(bf, af, data["FreeAcc_X"], method="gust")
        data["FreeAcc_Y_20hz"]=signal.filtfilt(bf, af, data["FreeAcc_Y"], method="gust")
        data["FreeAcc_Z_20hz"]=signal.filtfilt(bf, af, data["FreeAcc_Z"], method="gust")

        #calculate acc mag from filtered vectors
        data["Acc_magnitude_20hz"]=(data.FreeAcc_X_20hz**2 + data.FreeAcc_Y_20hz**2 + data.FreeAcc_Z_20hz**2)**(1/2)

        #save filtered dataset
        filename=ID+"-"+suf+"_"+"Acc_filter20hz_"+Phase+".csv"
        data_acc=data[["sample_number", "Acc_magnitude", "FreeAcc_X", "FreeAcc_Y", "FreeAcc_Z", "Acc_magnitude_20hz", "FreeAcc_X_20hz", "FreeAcc_Y_20hz", "FreeAcc_Z_20hz"]]

        os.chdir(OutputDirectory+Phase)
        data_acc.to_csv(filename)

        #Change back to data directory
        os.chdir(DataDirectory+Phase)

        print("filtered", ID, suf, Phase)