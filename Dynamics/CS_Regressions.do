** Housekeeping

*Import Data
clear all
cd "/Users/Maximilian/Documents/Princeton/RA Itskhoki/Itskhoki and Gaubert (2018)/Clean Code/Dynamics/Results/Data"
* import delimited regdata_357
import delimited CS_regdata_119

* Define auxiliary Variables
gen log_X = log(x)
gen log_D = log(d)

* Column 1
reg log_X top1 log_D

* Column 2
reg log_X top3 log_D

* Column 3
reg log_X top1 log_D t20share

* Column 4
reg log_X top3 log_D t20share
