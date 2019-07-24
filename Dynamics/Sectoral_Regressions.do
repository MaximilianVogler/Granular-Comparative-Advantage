** Housekeeping

*Import Data
clear all
cd "/Users/Maximilian/Documents/Princeton/RA Itskhoki/Itskhoki and Gaubert (2018)/Clean Code/Dynamics/Results/Data"
* import delimited regdata_357
import delimited sectoral_regdata_119

* Define auxiliary Variables
gen log_xt = log(x_t)
gen log_xt1 = log(x_t1)

eststo clear

** Table 1

* Column 1
eststo baseline: reg log_xt log_xt1

* Column 2
eststo time: reg log_xt log_xt1 i.year

* Column 3
* eststo sector: reg log_xt log_xt1 i.id

* Column 4
*eststo both: reg log_xt log_xt1 i.year i.id

esttab using "sectoral_119", replace star(* 0.10 ** 0.05 *** 0.01) booktabs b(a2) label

* Check FE (should give the same result as Column 4)
xtset id year

* Column 4
xtreg log_xt log_xt1, fe

* Column 5
xtreg log_xt log_xt1 i.year, fe
