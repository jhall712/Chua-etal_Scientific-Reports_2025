#obtain spectral power

%reset -f 

import os
import pandas as pd
import numpy as np

#define powersimp function

#create function to calculate power - loop over all infants still-face data, acc x, y, z, mag
from scipy import signal
from scipy.integrate import simps

def powersimp(ts, f_sampling, f_cutoff):

    # Define window length (4 seconds)
    sampling_rate=f_sampling
    win = 4 * sampling_rate

    # Compute psd
    data=ts
    freqs, psd = signal.welch(data, sampling_rate, nperseg=win)

    # Define delta lower and upper limits
    cutoff=f_cutoff
    low, high = 0, cutoff

    # Frequency resolution
    freq_res = freqs[1] - freqs[0]  # = 1 / 4 = 0.25

    # Find intersecting values in frequency vector
    idx_cutoff = np.logical_and(freqs >= low, freqs <= high)

    # Compute the absolute power by approximating the area under the curve
    power = simps(psd[idx_cutoff], dx=freq_res)

    # Relative power (expressed as a percentage of total power)
    total_power = simps(psd, dx=freq_res)
    rel_power = power / total_power

    return total_power, power, rel_power

#set directory to still-face segmented data directory
DataDirectory=r"XXXXX"
OutputDirectory=r"XXXXX"

os.chdir(DataDirectory)
Datafiles=os.listdir(DataDirectory)

row=0

#loop across all datafiles
for f in Datafiles:
    
    #skip invalid filenames
    if f.startswith("._"): continue
    elif f.endswith(".csv")==False: continue
    
    #obtain ID information - in first 4 characters of filename
    ID=f[0:4]
    
    #obtain suffix information
    #Specify file suffix information
    suffix= ["Wrist-Right",
         "Wrist-Left",
         "Ankle-Right",
         "Ankle-Left",
         "Torso"]
    
    for suf in suffix:
        if suf in f: break
    
        #if ID=="XXXX" and suf=="Torso":
    
        #Load data
        data=pd.read_csv(f, engine="python")

    for f_cutoff in [10, 20, 25, 30]:

        #loop across FreeAcc X, Y, Z to calculate total power, power within cutoff and relative power
        output_ID=pd.Series({"ID": int(ID), "f_cutoff": f_cutoff, "Sensor": suf})

        for axis in ["X", "Y", "Z"]:

            FreeAcc_axis="FreeAcc_"+axis
            ts=data[FreeAcc_axis]


            total_power, cutoff_power, cutoff_rel_power = powersimp(ts, 100, f_cutoff)

            output_cols={axis+"_total_power": total_power, 
                         axis+"_cutoff_power": cutoff_power, 
                         axis+"_cutoff_rel_power": cutoff_rel_power*100}

            if axis=="X": output_row=output_ID.append(pd.Series(output_cols))
            else: output_row=output_row.append(pd.Series(output_cols))

        output_row=pd.DataFrame(output_row).transpose()

        if row==0: df_output=output_row; row=1
        else: df_output=df_output.append(output_row)
            
    print(ID, suf, "- spectral power analysis completed")

#save dataframe with power information
df_output.reset_index(drop=True)
df_output

os.chdir(OutputDirectory)
df_output.to_csv("dataset-sensors-stillface-spectralpower.csv")