#obtain segmentation information

import os
import pandas as pd
import numpy as np

#Set directory
MetadataLocation=r"xxxxx"
OutputLocation=r"xxxxx"
DataLocation=r"xxxxx"

#Load SF phase info dataset
os.chdir(MetadataLocation)
SF_info=pd.read_csv("dataset-sensors-synchronisation-dec2022.csv")

SF_info["P_start"]=SF_info["P"]
SF_info["SF1_start"]=SF_info["SF1"]
SF_info["R1_start"]=SF_info["R1"]
SF_info["SF2_start"]=SF_info["SF2"]
SF_info["R2_start"]=SF_info["R2"]

#obtain SF phases start and end times
row_p1=7
row_p2=14
row_r1=11
row_r2=29
row_s1=34
row_r3=34

index=SF_info["P"][row_p1].index("-")
SF_info["P_start"].iloc[row_p1]=SF_info["P"][row_p1][:index]
print(SF_info["P_start"][row_p1])

index=SF_info["P"][row_p2].index("-")
SF_info["P_start"].iloc[row_p2]=SF_info["P"][row_p2][:index]
print(SF_info["P_start"][row_p2])

index=SF_info["R1"][row_r1].index("-")
SF_info["R1_start"].iloc[row_r1]=SF_info["R1"][row_r1][:index]
print(SF_info["R1_start"][row_r1])

index=SF_info["R1"][row_r2].index("-")
SF_info["R1_start"].iloc[row_r2]=SF_info["R1"][row_r2][:index]
print(SF_info["R1_start"][row_r2])

index=SF_info["SF2"][row_s1].index("-")
SF_info["SF2_start"].iloc[row_s1]=SF_info["SF2"][row_s1][:index]
print(SF_info["SF2_start"][row_s1])

index=SF_info["R2"][row_r3].index("-")
SF_info["R2_start"].iloc[row_r3]=SF_info["R2"][row_r3][:index]
print(SF_info["R2_start"][row_r3])


SF_info=SF_info.astype({"P_start": float, "SF1_start" : float, "R1_start":float, "SF2_start":float, "R2_start":float})

SF_info["P_end"]=SF_info["P_start"]+120
SF_info["SF1_end"]=SF_info["SF1_start"]+120
SF_info["R1_end"]=SF_info["R1_start"]+120
SF_info["SF2_end"]=SF_info["SF2_start"]+120
SF_info["R2_end"]=SF_info["R2_start"]+120

index=SF_info["P"][row_p1].index("-")
SF_info["P_end"].iloc[row_p1]=SF_info["P"][row_p1][index+2:]
print(SF_info["P_end"][row_p1])

index=SF_info["P"][row_p2].index("-")
SF_info["P_end"].iloc[row_p2]=SF_info["P"][row_p2][index+2:]
print(SF_info["P_end"][row_p2])

index=SF_info["R1"][row_r1].index("-")
SF_info["R1_end"].iloc[row_r1]=SF_info["R1"][row_r1][index+2:]
print(SF_info["R1_end"][row_r1])

index=SF_info["R1"][row_r2].index("-")
SF_info["R1_end"].iloc[row_r2]=SF_info["R1"][row_r2][index+2:]
print(SF_info["R1_end"][row_r2])

index=SF_info["SF2"][row_s1].index("-")
SF_info["SF2_end"].iloc[row_s1]=SF_info["SF2"][row_s1][index+2:]
print(SF_info["SF2_end"][row_s1])

index=SF_info["R2"][row_r3].index("-")
SF_info["R2_end"].iloc[row_r3]=SF_info["R2"][row_r3][index+2:]
print(SF_info["R2_end"][row_r3])


#Get number of samples for all infants

#Set directory
Folder=DataLocation

IDs=os.listdir(Folder)


print(IDs)
series=[]
ID_series=[]

for ID in IDs:
    os.chdir(Folder+ID)
    Datafiles=os.listdir(Folder+ID)
    for file in Datafiles:
        if file.startswith(ID) & file.endswith(".csv"):  
            print(file)
            df=pd.read_csv(file, skiprows=4, encoding="utf8", engine="python")
            samples=df.index.stop+1
            ID_series=ID_series + [ID]
            series=series + [samples]
            print(ID, samples)
            break
        
#Merge sample information with synchronisation dataset

s1=pd.Series(series, name='samples')
s2=pd.Series(ID_series, name='ID')

samples_id=pd.concat([s2, s1], axis=1)
samples_id["ID"]=samples_id["ID"].astype(int)

df=pd.merge(samples_id, SF_info, on="ID")
df=df.set_index(df["ID"])
df


#Calculate sensor start times for infants with only sensor synchronisation end times
sampling_rate=100

ID1=XXXX
ID2=XXXX

Sensor_end=df["Sensor_end"][ID1]
samples=df["samples"][ID1]

Sensor_start=Sensor_end-samples/sampling_rate

df["Sensor_start"].loc[ID1]=Sensor_start

Sensor_end=df["Sensor_end"][ID2]
samples=df["samples"][ID2]

Sensor_start=Sensor_end-samples/sampling_rate

df["Sensor_start"].loc[ID2]=Sensor_start

df

#Calculate sensor data segments in sample time
sampling_rate=100

#Check data type
df=df.astype({"P_start": float, "SF1_start" : float, "R1_start":float, "SF2_start":float, "R2_start":float})
df=df.astype({"P_end": float, "SF1_end" : float, "R1_end":float, "SF2_end":float, "R2_end":float})

#Convert video times to sensor time
df["P_start_seg"]=(df["P_start"]-df["Sensor_start"])*sampling_rate
df["P_end_seg"]=(df["P_end"]-df["Sensor_start"])*sampling_rate

df["SF1_start_seg"]=(df["SF1_start"]-df["Sensor_start"])*sampling_rate
df["SF1_end_seg"]=(df["SF1_end"]-df["Sensor_start"])*sampling_rate

df["R1_start_seg"]=(df["R1_start"]-df["Sensor_start"])*sampling_rate
df["R1_end_seg"]=(df["R1_end"]-df["Sensor_start"])*sampling_rate

df["SF2_start_seg"]=(df["SF2_start"]-df["Sensor_start"])*sampling_rate
df["SF2_end_seg"]=(df["SF2_end"]-df["Sensor_start"])*sampling_rate

df["R2_start_seg"]=(df["R2_start"]-df["Sensor_start"])*sampling_rate
df["R2_end_seg"]=(df["R2_end"]-df["Sensor_start"])*sampling_rate

df_seg=df[["samples", "P_start_seg", "P_end_seg", "SF1_start_seg", "SF1_end_seg", "R1_start_seg","R1_end_seg","SF2_start_seg", "SF2_end_seg", "R2_start_seg","R2_end_seg"]]

df_seg=round(df_seg,0)

df_seg=df_seg.sort_index()

#save segmentation info
os.chdir(MetadataLocation)
df_seg.to_csv("dataset-sensors-segmentation-dec2022.csv")