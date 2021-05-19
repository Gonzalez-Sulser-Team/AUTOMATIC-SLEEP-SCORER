# AUTOMATIC-SLEEP-SCORER
AUTOMATIC SLEEP SCORING PIPELINE - Dr. Ingrid Buller &amp;  Dr. Alejandro Bassi  for  @Gonzalez-Sulser-Team
### AUTOMATIC EEG AND SLEEP SCORING PIPELINE INSTRUCTIONS

**_Dr. Ingrid Buller - ingrid.buller@ed.ac.uk_**

**_Dr. Alejandro Bassi abassi@dcc.uchile.cl - Dr. Adrian Ocampo-Garces aocampo@med.uchile.cl_**

**_@Gonzalez-Sulser-Team (in collaboration with Sleep & Chronobiology Lab - ICBM - Faculty of Medicine - Universidad de Chile)_**

**_CDBS - SIDB - University of Edinburgh // ICBM - Faculty of Medicine - Universidad de Chile_**

### **STEP ONE: OPEN .DAT FILE AND CHOOSE CHANNELS FOR SLEEP ANALYSIS (FOR WIRELESS RECORDING OBTAINED THROUGH TAINI SYSTEM) - OPTIONS:**

**_A) OPEN FULL .DAT FILE AND SELECT ONE EEG AND ONE EMG CHANNEL:_**
Data parsing scripts provided by the TAINITEC team can be found here (to open all 16 channels) 
https://bitbucket.org/tainitec/data_parsing_samples/src/master/

a.1) Set start and end samples number to obtain the full 24 hr-day (starting at 07:00:00 to 06:59:59 on the following day)

a.2) Choose an EEG channel (4 or 13 suggested) and an EMG channel (2 or 15) for the sleep scoring analysis.

**_B) CHECK AND CHOOSE BETWEEN 2 EEG AND 2 EMG CHANNELS, AND SAVE INTO .CSV FILE (**MATLAB_FILES folder, Script: Taini_CheckChoose.m**)_**

b.1) To run a set of animals create a .csv file using the format attached as **rat_day_index_CDKL5.csv** (Note BL3 in this case are split in
two files that need to be merged, the script is set to do this automatically). 
By default columns M and O are set to "0" and can be automatically filled by the script with the selected eeg and emg channel numbers.

b.2) Change names for "rat", "day", "line", "pn", "p" variables (Note source .dat files are set as p '\' line '\' ratname),
with ratname = [line '_' rat]. (Example for file and Folder format for data: source_folder/CDKL5/CDKL5_1959/)

b.3) Run the script. 

b.4) 2 options of EEG channels will appear along with a prompt window asking to choose (insert 1 or 2).

b.5) 2 options of EMG channels will appear along with a prompt window asking to choose (insert 1 or 2).

b.6) A .cvs file will be written containing the chosen EEG and EMG channels with the format [ratname '_' day '_Channels.csv'], into a folder
that can be specified in line 178 of the script.

### **OPTIONAL STEP: SWDs AUTOMATIC DETECTION (MATLAB ONLY)**

1) A prompt window will ask if you would like to run SWDs auto-detection scripts (**TAINI_autodetect.m** and **cepstral_analysis.m** plus additional
scripts included in **MATLAB_FILES folder**)

2) After SWD analysis is done, a prompt window will ask you if you want to check the events. **_ALWAYS ANSWER YES AND CHECK IF THE EVENTS SHOWN 
ARE SWD IN DEED._**

3) If SWDs are correctly detected, choose "YES" to save results, column L on rat_day_index_CDKL5.csv will be automatically filled with 1 if SWDs are found.

### **STEP TWO: RUN AUTOMATIC SLEEP SCORER ON R STUDIO (**R_SleepAutoscorer_FILES folder**, script: Sleep_autoscore_Beta.R)**

**_Original version by Dr. Alejandro Bassi and Dr. Javier Díaz in https://doi.org/10.3389/fncel.2017.00302  - Sleep & Chronobiology Lab - ICBM - Faculty of Medicine - Universidad de Chile_**

**_Ported to R by Dr. Alejandro Bassi (abassi@dcc.uchile.cl)_**

**_Adapted by Dr. Ingrid Buller (ingrid.buller@ed.ac.uk)_**

**_@Gonzalez-Sulser-Team CDBS - SIDB - University of Edinburgh 2021_**

_A) Paste all the files inside the **R_SleepAutoscorer_FILES folder (sleep_score_library.R , mgg.R , Sleep_autoscore_Beta.R)** into the source folder where you will keep your data (i.e "C:\Users\yourname\Documents\POSTDOC\ANALISYS\EEG RESULTS\Sleep_Results")_

_B) ONCE WHEN OPENING R-STUDIO FOR NEW SESSION:_

b.1) When opening for the first time: Install all required packages from CRAN (**dplyr psych multitaper rgl ggplot2**)

b.2) Set corresponding directories: line 21 SET SOURCE DIR FOR CODE AND LIBRARY, 42 (source dir for eeg and emg channels)

b.3) Type in console:

setwd("/Users/yourname/Documents/POSTDOC/ANALISYS/EEG RESULTS/Sleep_Results") ###Example source dir

source("mgg.R")

source("sleep_score_library.R")

library(rgl)

source('~/POSTDOC/ANALISYS/EEG RESULTS/Sleep_Results/Sleep_autoscore_Beta.R') #### example script source folder

**_AND TYPE VARIABLES TO CHANGE FOR EACH RAT/DAY:_** 

**rat_line<<-("CDKL5") ## CHANGE AND PASTE TO CONSOLE**

**animal_id<<-("1959") ## CHANGE AND PASTE TO CONSOLE**

**animal_day<<-("BL1") ## CHANGE AND PASTE TO CONSOLE**

**seizures<<-0 ##RECORD WITH SEIZURES? 1=YES 0=NO**

(Example for file and Folder format for data: source_folder/CDKL5/CDKL5_1959/)

_C) RUN THE MAIN FUNCTION TO OBTAIN AUTOMATIC SLEEP SCORING AND POWER SPECTRUM BY STATE_
**sleep.autoscore()**

_OR RUN THE FOLLOWING FUNCTIONS JUST TO OBTAIN AUTOMATIC SLEEP SCORING WITHOUT POWER SPECTRUM BY STATE_

**rat.get()**

**rat.correct()**

_D) OUTPUTS:_ 

d.1) Line 178: .csv with sleep scoring only (0=WAKE, 1=NON REM, 2=REM, EPOCHS 5 SECS)

d.2) Line 186: .csv with sleep scoring and swd for coherence analysis (0=WAKE, 1=NON REM, 2=REM, 4=SWD, EPOCHS 5 SECS).

d.3) Line 272: .csv with power spectrum by state - Column names: s_0=Wakefulness, s_1=Non-REM, s_2=REM (if SWD, s_4=SWD). 

d.4) 3D graph showing the initial distribution of the epochs in the axes: Theta power, EEG (-theta) power, and EMG power.
(Green = Wakefulness, Blue = Non-REM Sleep, Red = REM sleep, Black = epochs still to be classified by context)

d.5) EEG spectrogram in 12 consecutive rows of 2 hours for the range 0.4-20 Hz (Blue to Yellow = low to high power)

d.6) Distribution of states by hour of recording in epocs (one hour = 720 epochs of 5 secs). Green = Wakefulness, Blue = Non-REM Sleep, Red = REM sleep. 
