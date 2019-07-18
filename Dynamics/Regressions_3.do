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
eststo CS_2005: reg log_x top3 log_d if year==2005

* Column 3
eststo Panel: reg log_x top3 log_d i.year , vce(cluster id)

* Column 3'
eststo Panel_prime: reg log_x top3 i.year , vce(cluster id)

* Column 5
xtset id year
eststo Dynamics_1: xtreg log_x top3 i.year, fe 

* Column 6
eststo Dynamics_2: reg d.log_x d.top3


keep if year == 1997 | year == 2008
sort id year
gen Delta_log_x = log_x[_n+1]-log_x[_n]
gen Delta_top_3 = top3[_n+1]-top3[_n]
keep if year == 1997

* Column 8
eststo Dynamics_3: reg Delta_log_x Delta_top_3

esttab using "table_1_big", replace star(* 0.10 ** 0.05 *** 0.01) booktabs b(a2) label


** Table 2

eststo clear

* Column 1
eststo OLS_1: reg Delta_log_x log_x log_d 

* Column 2
eststo OLS_2: reg Delta_log_x top3 log_d

* Column 3
eststo OLS_3: reg Delta_log_x log_x top3 log_d

* Column 4
eststo IV_1: ivregress 2sls Delta_log_x log_d (log_x = top3)

* Column 5
eststo IV_2: ivregress 2sls Delta_log_x (log_x = top3)

esttab using "table_2_big", replace star(* 0.10 ** 0.05 *** 0.01) booktabs b(a2) label
