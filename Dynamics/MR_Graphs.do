***** Appendix Tables *****


** Housekeeping

*Import Data
clear all
cd "/Users/Maximilian/Documents/Princeton/RA Itskhoki/Itskhoki and Gaubert (2018)/Clean Code/Dynamics/Results/Data"
import delimited regdata_9996_extra

* Define auxiliary Variables
gen log_x = log(x)
gen log_d = log(d)

eststo clear

keep if year == 0 | year == 50
sort id year
gen Delta_log_x = log_x[_n+1]-log_x[_n]
gen Delta_top_1 = top1[_n+1]-top1[_n]
gen Delta_top_3 = top3[_n+1]-top3[_n]
keep if year == 0

* Column 2 for Top 1
eststo OLS_2_1: reg Delta_log_x top1 log_d

* Column 2 for Top 3
eststo OLS_2_3: reg Delta_log_x top3 log_d 

* Column 3 for Top 1
eststo OLS_3_1: reg Delta_log_x log_x top1 log_d

* Column 3 for Top 3
eststo OLS_3_3: reg Delta_log_x log_x top3 log_d


esttab using "table_Y50", replace star(* 0.10 ** 0.05 *** 0.01) booktabs b(a2) label
