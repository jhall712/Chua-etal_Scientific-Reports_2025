import os
import pandas as pd
import numpy as np

#Set directory
MetadataLocation=r"XXXXX"
OutputLocation=r"XXXXX"
DataLocation=r"XXXXX"

#Load segmentation information
os.chdir(MetadataLocation)
segments=pd.read_csv("dataset-sensors-segmentation-dec2022.csv")

#Obtain list of IDs to segment
IDs=segments.ID.astype(str)

#Specify file suffix information
suffix= ["-Wrist-Right",
         "-Wrist-Left",
         "-Ankle-Right",
         "-Ankle-Left",
         "-Torso"]

#Specify phases to segment
phases=["P", "SF1", "R1", "SF2", "R2"]


segments=segments.set_index(segments["ID"])
segments

for ID in IDs:
    
    #set directory and obtain list of filenames
    os.chdir(DataLocation+ID)
    Datafiles=os.listdir(DataLocation+ID)
    
#    if ID=="8208":
    for f in Datafiles:
        if f.startswith("._"): continue
        elif f.endswith(".csv")==False: continue

        for phase in phases:
            skip_phase="No"
            col_start=phase+"_start_seg"
            col_end=phase+"_end_seg"

            time_start=segments[col_start].loc[int(ID)]
            time_end=segments[col_end].loc[int(ID)]

            if pd.isna(time_start): skip_phase="Yes"

            if skip_phase=="Yes": 
                print("skip", ID, phase, f)
                continue


            for suf in suffix:
                file_suf=ID+suf

                #read filename based on phase and suffix
                if f.startswith(file_suf):

                    #segment data
#                        if phase=="P" and suf=="-Torso":
#                    print("segment", ID, suf, phase)

                    df=pd.read_csv(f, skiprows=4, engine="python")

                    dfslc=df.iloc[int(time_start): int(time_end)]

                    #obtain acceleration magnitude
                    acc=(dfslc.FreeAcc_X**2 + dfslc.FreeAcc_Y**2 + dfslc.FreeAcc_Z**2)**(1/2)

                    #obtain gyroscope magnitude
                    gyr=(dfslc.Gyr_X**2 + dfslc.Gyr_Y**2 + dfslc.Gyr_Z**2)**(1/2)

                    #obtain segmented dataframe
                    s1=pd.Series(acc, name='Acc_magnitude')
                    s2=pd.Series(gyr, name='AngVel_magnitude')

                    s3=dfslc[["FreeAcc_X", "FreeAcc_Y", "FreeAcc_Z", "Gyr_X", "Gyr_Y", "Gyr_Z"]]

                    df_seg=pd.concat([s1, s2, s3], axis=1)

                    df_seg=df_seg.reset_index()
                    df_seg=df_seg.rename(columns={"index": "sample_number"})

                    #save dataset in corresponding directory
                    OutputDir=OutputLocation+"/"+phase
                    os.chdir(OutputDir)

                    outputname=ID + suf + "_" + phase +".csv"
                    df_seg.to_csv(outputname)

                    #change back to data directory
                    os.chdir(DataLocation+ID)
                    
    print(ID, " - segmenting completed")
    
df_seg.rename(columns={"index": "sample_number"})


#Segment entire still-face

#Obtain list of IDs to segment
IDs=segments.ID.astype(str)

for ID in IDs:
    
    #set directory and obtain list of filenames
    os.chdir(DataLocation+ID)
    Datafiles=os.listdir(DataLocation+ID)
    
#    if ID=="XXXX":
    for f in Datafiles:
        if f.startswith("._"): continue
        elif f.endswith(".csv")==False: continue

        for phase in ["Still-face"]:
            skip_phase="No"
            col_start="start"
            col_end="end"

            time_start=segments[col_start].loc[int(ID)]
            time_end=segments[col_end].loc[int(ID)]

            if pd.isna(time_start): skip_phase="Yes"

            if skip_phase=="Yes": 
                print("skip", ID, phase, f)
                continue


            for suf in suffix:
                file_suf=ID+suf

                #read filename based on phase and suffix
                if f.startswith(file_suf):

                    #segment data
#                        if phase=="P" and suf=="-Torso":
#                    print("segment", ID, suf, phase)

                    df=pd.read_csv(f, skiprows=4, engine="python")

                    dfslc=df.iloc[int(time_start): int(time_end)]

                    #obtain acceleration magnitude
                    acc=(dfslc.FreeAcc_X**2 + dfslc.FreeAcc_Y**2 + dfslc.FreeAcc_Z**2)**(1/2)

                    #obtain gyroscope magnitude
                    gyr=(dfslc.Gyr_X**2 + dfslc.Gyr_Y**2 + dfslc.Gyr_Z**2)**(1/2)

                    #obtain segmented dataframe
                    s1=pd.Series(acc, name='Acc_magnitude')
                    s2=pd.Series(gyr, name='AngVel_magnitude')

                    s3=dfslc[["FreeAcc_X", "FreeAcc_Y", "FreeAcc_Z", "Gyr_X", "Gyr_Y", "Gyr_Z"]]

                    df_seg=pd.concat([s1, s2, s3], axis=1)

                    df_seg=df_seg.reset_index()
                    df_seg=df_seg.rename(columns={"index": "sample_number"})

                    #save dataset in corresponding directory
                    OutputDir=OutputLocation+"/"+phase
                    os.chdir(OutputDir)

                    outputname=ID + suf + "_" + phase +".csv"
                    df_seg.to_csv(outputname)

                    #change back to data directory
                    os.chdir(DataLocation+ID)
                    
    print(ID, " - segmenting completed")