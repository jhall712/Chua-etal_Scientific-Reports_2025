#Run IMPE code across P, SF1, R1, SF2, R2 for selected parameters: m=4, tlag=1, scales 1-50
%reset -f 

import os
import pandas as pd
import numpy as np

SourceDirectory=r"XXXXX"
os.chdir(SourceDirectory)
from IMPE import pe, CoarseGrain, impe


#Specify file suffix information
suffix= ["Wrist-Right",
         "Wrist-Left",
         "Ankle-Right",
         "Ankle-Left",
         "Torso"]

#Set parameters
m=4
tlag=1
scales=list(range(1,51))

#obtain list of IDs 
#set metadata location
MetadataLocation=r"XXXXX"
os.chdir(MetadataLocation)
#metadata=pd.read_csv("dataset-sensors-segmentation.csv")
metadata=pd.read_csv("dataset-sensors-segmentation-dec2022.csv")
IDs=metadata.ID.tolist()
print(IDs)
#exclude=[XXXX, XXXX, XXXX, XXXX, XXXX, XXXX]
#for element in IDs:
#    if element in exclude:
#        print(element)
#        IDs.remove(element)
#IDs.remove(XXXX)
        
phase="R2"

#set data directory
DataDirectory=r"XXXXX"
FilteredData=DataDirectory + phase

filtered="yes"

files=os.listdir(DataDirectory + phase)

for filtered in ["yes"]:
    df=[]
    row=0
    os.chdir(FilteredData)
    for ID in IDs:
        for suf in suffix:
            filename=str(ID)+"-"+suf+"_Acc_filter20hz_"+phase+".csv"
            
            if filename not in files: continue

            data=pd.read_csv(filename, engine="python")

            if filtered=="yes":
                TS=data.Acc_magnitude_20hz

            elif filtered=="no":
                TS=data.Acc_magnitude

            TS=TS.tolist()

            IMPE_all=[np.nan for x in range(1, len(scales)+1+3)]

            IMPE_all[0]=ID
            IMPE_all[1]=phase
            IMPE_all[2]=suf

            #put into IMPE code with parameters embedding dimension, lag

            index=0
            for scale in scales:
                imPE=impe(TS, m, tlag, scale)
                IMPE_all[index+3]=imPE
                index=index+1

            if row==0:
                df=[IMPE_all]
                row=1
            else:
                df=df+[IMPE_all]

            print(ID, phase, "filtered", filtered, suf, "m=", m, "complete")

    #save IMPE output
    colnames=["ID", "Phase", "Sensor"] + scales
    IMPE_output=pd.DataFrame(df, columns=colnames)

    #file suffix
    filesuf="m"+str(m)+"lag"+str(tlag)+"scale50"

    if filtered=="yes": 
        output_filename="dataset-IMPE-" + "filtered" + "_" + phase + "_" + filesuf + ".csv"
    if filtered=="no": 
        output_filename="dataset-IMPE-" + "unfiltered" + "_" + phase + "_" + filesuf + ".csv"

    OutputDirectory=r"XXXXX"
    os.chdir(OutputDirectory)
    IMPE_output.to_csv(output_filename)