# __Magnetophosphene Perception in the ELF-MF range (0-300Hz)__
Lawson Health Research Institute - Imaging - Human Threshold Research Group
Dr. S. Villard 2020


This project look at the frequency response to Extremely Low Frequency Magnetic Field (ELF-MF) exposure.
The perception of flickering lights when closing your eyes in the presence of ELF-MF has fisrt been described by d'[Arsonval](https://en.wikipedia.org/wiki/Jacques-Arsène_d%27Arsonval) in 1896:  
 > d’Arsonval, A. (1896). Dispositifs pour la mesure des courants alternatifs de toutes fréquences. Compt. Rend. Soc. Biol, 3, 450--451.
 
It has been used since has the principal biomarker of the impact of ELF-MF on the brain. Estbalishment of threshold for magnetophosphene perception is at the core of international guidelines and standards such as the [ICNIRP](https://www.icnirp.org/cms/upload/publications/ICNIRPlfgaps2020.pdf) or [IEEE](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8910342).

In this project we investigated the perception threshold of healthy participants from 5 to 300Hz.

__We are trying to model the frequency response of the magnetophosphene perception__

## 1. Usage

 - main notebook: ```2020_04_FrequencyResponse.ipynb```
 
## 2. Analysis

 - compute Perception Rate for each frequency
 - extract Frequencies with > 80% of perception to build our model
 - visualization of raw data of perception threshold
 - regression model: We can choose to build the model over the perception express in Flux density and then convert to dB/dt or the opposite. Since our hypothesis is based on the induced electric field at the retinal level we chose to build our model on the dB/dt threshold values and then represent these values in terms of Flux density.