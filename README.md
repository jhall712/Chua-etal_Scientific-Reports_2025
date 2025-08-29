[README.md](https://github.com/user-attachments/files/22046557/README.md)

A multi-level analysis of motor and behavioural dynamics in 9-month-old preterm and term-born infants during changing emotional and interactive contexts


** Project overview **

This repository is linked to a project assessing the entropy of infant motor acceleration during the still-face paradigm. It contains the analysis dataset, motor data preprocessing scripts, and the scripts for statistical analyses. The data is based on a subsample of the Theirworld Edinburgh Birth Cohort (TEBC).

The still-face paradigm contains five phases (play, still-face 1, reunion 1, still-face 2, reunion 2)
Motor accelerometer data was collected from five sensors on the wrists (right, left), legs (right, left) and torso. Permutation entropy at 50 timescales (ie., Multiscale permutation entropy) was calculated from accelerometer data. Frequency-specific permutation entropy was obtained by summing permutation entropy within each frequency bands (gamma, beta, alpha, theta, delta).

Chua, Y., Jiménez-Sánchez, L., Ledsham, V. et al. A multi-level analysis of motor and behavioural dynamics in 9-month-old preterm and term-born infants during changing emotional and interactive contexts. Sci Rep 15, 952 (2025). https://doi.org/10.1038/s41598-024-83194-w

** Files **

Analysis dataset

TEBC-dataset-IMPE-filtered-m4lag1scale50-all-phases-with-demog.csv: This de-identified dataset contains the processed motor entropy variables, alongside other demographic information from the TEBC.

Preprocessing scripts folder

This folder contains Python scripts used for pre-processing raw tri-axial accelerometer data

1. SegmentationInfo.py: This script uses a file containing sensor and video synchronisation information to convert timestamps from the Still-face paradigm to timestamps in sensor time, so that sensor data periods can be obtained corresponding to the still-face paradigm.

2. SegmentSensorData.py: Accelerometer files from each sensor is segmented into five periods corresponding to each phase of the still-phase paradigm.

3. FilterSensorData.py: Applies low pass 20Hz filter

4. ObtainSpectralPower.py: Calculates the spectral power 

5. IMPE.py: This file defines the functions to obtain multiscale permutation entropy using the Improved Multiscale Permutation Entropy (IMPE) algorithm, including the permutation entropy algorithm, a coarsegraining function, and the IMPE function which calculates permutation entropy after coarsegraining

6. IMPE-allPhases-loop.py: Script to loop through segmented data files to calculate IMPE (Azami & Escudero, 2016) for filtered and unfiltered data according to input parameters m, tlag and scale

7. IMPE-sensitivityAnalyses-loop.py: Calculates IMPE according to a different embedding dimension to compare sensitivity of IMPE to parameter selection.

Azami, H. & Escudero, J. Improved multiscale permutation entropy for biomedical signal analysis: Interpretation and application to electroencephalogram recordings. Biomed. Signal Process. Control 23, 28–41 (2016).


Analysis scripts

1. SampleDescription-entropy.R: runs descriptive analyses on the analysis dataset
 
2. SampleDescription-entropy-additional.R: runs descriptive analyses in relation to infants sitting ability

3. Analysis-Entropy-GrowthCurves.R: runs linear mixed effect models to analyse effect of Preterm birth, Sensor and Scale on permutation entropy

4. Analysis-Entropy-SF-paradigm.R: runs linear mixed effect models to analyse effect of Preterm birth, Sensor and Phase on permutation entropy in gamma, beta, alpha, theta, delta frequency bands

5. entropy-sensitivity-analyses.R: re-runs statistical analyses excluding infants who are not able to sit independently

** Contributors **
Yu Wei Chua, Lorena Jiménez Sánchez, Victoria Ledsham, Sinéad O’Carroll, Ralf F. A. Cox, Ivan Andonovic, Christos Tachtatzis, James P. Boardman, Sue Fletcher-Watson, Philip Rowe and Jonathan Delafield-Butt

** License **
CC-BY 4.0
http://creativecommons.org/licenses/by/4.0/
