#Run IMPE code across SF1 unfiltered and filtered data
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
M=[5]
tlag=1
scales=list(range(21,51))

#obtain list of IDs 
#set metadata location
MetadataLocation=r"XXXXX"
os.chdir(MetadataLocation)
metadata=pd.read_csv("dataset-sensors-segmentation.csv")
IDs=metadata.ID.tolist()
exclude=[XXXX, XXXX, XXXX, XXXX]
for element in IDs:
    if element in exclude:
        print(element)
        IDs.remove(element)

phase="SF1"

#set data directory
FilteredData=r"XXXXX"

for m in M:
    for filtered in ["yes", "no"]:
        df=[]
        row=0
        os.chdir(FilteredData)
        for ID in IDs:
            for suf in suffix:
                filename=str(ID)+"-"+suf+"_Acc_filter20hz_"+phase+".csv"

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

                print(ID, "SF1 filtered", filtered, suf, "m=", m, "complete")

        #save IMPE output
        colnames=["ID", "Phase", "Sensor"] + scales
        IMPE_output=pd.DataFrame(df, columns=colnames)

        #file suffix
        filesuf="m"+str(m)+"lag"+str(tlag)+"scale50"

        if filtered=="yes": 
            output_filename="dataset-IMPE-" + "filtered" + "_SF1_" + filesuf + ".csv"
        if filtered=="no": 
            output_filename="dataset-IMPE-" + "unfiltered" + "_SF1_"+ filesuf + ".csv"

        OutputDirectory=r"XXXXX"
        os.chdir(OutputDirectory)
        IMPE_output.to_csv(output_filename)