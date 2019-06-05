***** Appendix Tables *****


** Housekeeping

*Import Data
clear all
cd "/Users/Maximilian/Documents/Princeton/RA Itskhoki/Itskhoki and Gaubert (2018)/Clean Code/Dynamics/Results/Data"
* import delimited regdata_357
import delimited regdata_9996

* Define auxiliary Variables
gen log_x = log(x)
gen log_d = log(d)

eststo clear

** Table 1

* Column 1
eststo CS_2005: reg log_x top1 log_d if year==2005, vce(cluster id)

* Column 3
xtset id year
eststo Panel: xtreg log_x top1 log_d i.year , fe vce(cluster id)

* Column 5
* xtreg log_x top1 i.year i.id, fe vce(cluster id)

* Column 6
eststo Dynamics: xtreg d.log_x d.top1, vce(cluster id)

esttab using "table_3", replace star(* 0.10 ** 0.05 *** 0.01) booktabs b(a2) label
** Table 2

eststo clear

* Column 1
keep if year == 1997 | year == 2008
sort id year
gen Delta_log_x = log_x[_n+1]-log_x[_n]
keep if year == 1997

eststo OLS_1: reg Delta_log_x log_x log_d 

* Column 2
eststo OLS_2: reg Delta_log_x top1 log_d

* Column 3
eststo OLS_3: reg Delta_log_x log_x top1 log_d

* Column 4
eststo IV_1: ivregress 2sls Delta_log_x log_d (log_x = top1)

* Column 5
eststo IV_2: ivregress 2sls Delta_log_x (log_x = top1)

esttab using "table_4", replace star(* 0.10 ** 0.05 *** 0.01) booktabs b(a2) label
