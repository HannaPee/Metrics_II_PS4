//// Metrics PS 4 /// 

** housekeeping
clear all                   // remove anything old stored
set more off, permanently   // tell Stata not to pause
set linesize 255            // set line length for the log file
version                     // check the version of the command interpreter

* Set working directory to the current repo folder
cd "C:\Users\42610\OneDrive - Handelshögskolan i Stockholm\Documents\Metrics_II_PS4"
global wd "`c(pwd)'"

* Create folders if they do not exist
cap mkdir figures
cap mkdir output
cap mkdir logs

**** Q1 ****

*** 5 *** 

clear all 
set more off 
set obs 100 

* Create village ID
gen village_id = _n

* Expand each village to 10 individuals
expand 10

* Create individual ID within village
bysort village_id: gen individual_id = _n

** Generate treatment variable on village level***

set seed 12345

* Keep one row per village
preserve
collapse (count) n=individual_id, by(village_id)

* Randomize villages
gen u = runiform()
sort u

* Assign treatment (e.g., 50 treatment, 50 control)
gen D_j = (_n <= 50)
drop u n 

* Randomize mu * 
gen mu_j = rnormal(2,1)

* Save assignment
tempfile village_assign
save `village_assign'

restore

* Merge back to individual-level data
merge m:1 village_id using `village_assign'
drop _merge

* Create individual-level treatment within villages
gen u_i = runiform()
bysort village_id (u_i): gen W_ij = (_n <= 5)
drop u_i 

** Generate potential outcomes ** 

gen Y_00 = rnormal()
gen Y_10 = rnormal(1,1)
gen Y_01 = rnormal(mu_j,1)


*** 6 *** 
preserve 

* Create results file
tempname h1 h2
tempfile spec1 spec2
postfile `h1' phi b_d ll_d ul_d b_w ll_w ul_w using `spec1', replace
postfile `h2' phi b_d ll_d ul_d b_w ll_w ul_w b_dw ll_dw ul_dw using `spec2', replace

forvalues i = -10/10 {
    local phi = `i'/10

    capture drop X_ij Y_11 Y

    gen X_ij = rnormal(`phi', sqrt(.2))
    gen Y_11 = Y_10 + Y_01 + X_ij
    gen Y = (1-D_j)*(1-W_ij)*Y_00 + (1-D_j)*W_ij*Y_01 + D_j*(1-W_ij)*Y_10 + D_j*W_ij*Y_11

    * Run regression 1
    regress Y D_j W_ij, vce(cluster village_id)

    local tcrit = invttail(e(df_r), .025)

    local bd   = _b[D_j]
    local bdl  = _b[D_j] - `tcrit'*_se[D_j]
    local bdu  = _b[D_j] + `tcrit'*_se[D_j]

    local bw   = _b[W_ij]
    local bwl  = _b[W_ij] - `tcrit'*_se[W_ij]
    local bwu  = _b[W_ij] + `tcrit'*_se[W_ij]

    post `h1' (`phi') (`bd') (`bdl') (`bdu') (`bw') (`bwl') (`bwu')

    * Run regression 2
    regress Y D_j W_ij c.D_j#c.W_ij, vce(cluster village_id)

    local tcrit = invttail(e(df_r), .025)

    local bd   = _b[D_j]
    local bdl  = _b[D_j] - `tcrit'*_se[D_j]
    local bdu  = _b[D_j] + `tcrit'*_se[D_j]

    local bw   = _b[W_ij]
    local bwl  = _b[W_ij] - `tcrit'*_se[W_ij]
    local bwu  = _b[W_ij] + `tcrit'*_se[W_ij]

    local bdw  = _b[c.D_j#c.W_ij]
    local bdwl = _b[c.D_j#c.W_ij] - `tcrit'*_se[c.D_j#c.W_ij]
    local bdwu = _b[c.D_j#c.W_ij] + `tcrit'*_se[c.D_j#c.W_ij]

    post `h2' (`phi') (`bd') (`bdl') (`bdu') (`bw') (`bwl') (`bwu') (`bdw') (`bdwl') (`bdwu')
}

postclose `h1'
postclose `h2'


* First graph* 

*store true treatment effects* 
local tau_d = 1
local tau_w = 2

use `spec1', clear
sort phi

twoway ///
    (rcap ll_d ul_d phi, sort) ///
    (connected b_d phi, sort msymbol(O)) ///
    (rcap ll_w ul_w phi, sort) ///
    (connected b_w phi, sort msymbol(D)) ///
    , ///
    yline(`tau_d', lpattern(dash)) ///
    yline(`tau_w', lpattern(dot)) ///
    xtitle("phi") ///
    ytitle("Coefficient / true value") ///
    legend(order(2 "b(Dj)" 4 "b(Wij)" 5 "true tau_d" 6 "true tau_w"))
	
graph export "figures/reg1.pdf", replace 
	
* Second graph* 

use `spec2', clear

twoway ///
    (rcap ll_d ul_d phi, sort) ///
    (connected b_d phi, sort msymbol(O)) ///
    (rcap ll_w ul_w phi, sort) ///
    (connected b_w phi, sort msymbol(D)) ///
    (function y=`tau_d', range(-1 1) lpattern(dash)) ///
    (function y=`tau_w', range(-1 1) lpattern(dot)), ///
    xtitle("phi") ytitle("Coefficient / true value") ///
    legend(order(2 "b(Dj)" 4 "b(Wij)" 5 "true tau_d" 6 "true tau_w"))


graph export "figures/reg2.pdf", replace 

restore 

*** 7 *** 
* Generate variables* 
local phi = 0 
gen X_ij = rnormal(`phi', sqrt(.2))
gen Y_11 = Y_10 + Y_01 + X_ij
gen Y = (1-D_j)*(1-W_ij)*Y_00 + (1-D_j)*W_ij*Y_01 + D_j*(1-W_ij)*Y_10 + D_j*W_ij*Y_11

* Regress Y on Wij with village fixed effects* 

areg Y W_ij, absorb(village_id) vce(cluster village_id)

* Regress Y on W_ij and D_j * 

reg Y W_ij D_j, vce(cluster village_id)


*** 8 *** 

* Run regression
regress Y i.village_id##c.W_ij, vce(cluster village_id)

* Create dataset with one observation per village

gen tauhat = .

forvalues j = 1/100 {
    quietly {
        * average effect
        local base = _b[W_ij]

        * interaction term 
        capture local dev = _b[`j'.village_id#c.W_ij]

        if _rc != 0 local dev = 0

        replace tauhat = `base' + `dev' if village_id == `j'
    }
}

* Create a scatterplot * 

preserve
collapse (mean) mu_j, by(village_id)

twoway scatter tauhat mu_j, ///
    xtitle("μ_j") ///
    ytitle("Estimated τ_j") 

graph export "figures/tau_j.pdf", replace

 
restore 

*** 9 ***
*Calculate tau_eb** 

local sigma2_Y = 1
local sigma2_mu = 1

gen tau_eb = (4*`sigma2_Y'/10)/(`sigma2_mu' + 4*`sigma2_Y'/10)*mu_j + ///
             (`sigma2_mu')/(`sigma2_mu' + 4*`sigma2_Y'/10)*tauhat
			 
twoway ///
    (scatter tauhat mu_j, mcolor(gs12)) ///
    (scatter tau_eb mu_j, mcolor(navy)) ///
    (pcarrow tauhat mu_j tau_eb mu_j, lcolor(gs8)) ///
    (function y=x, range(-1 4) lpattern(dash)), ///
    xtitle("True μ_j") ///
    ytitle("Estimate of τ_j") ///
    legend(order(1 "Noisy" 2 "EB shrunk" 4 "45° line"))
			 

graph export "figures/tau_eb.pdf", replace 



*** Question 2 ****

**** Simulate data ****

clear all                 
set more off 

set seed 240505

* Dimensions
local N = 10
local T = 200
local J = 20

* 1) Draw worker effects: alpha_i ~ N(0,1)
preserve
clear
set obs `N'
gen worker = _n
gen alpha  = rnormal(0,1)
tempfile workers
save `workers'
restore

* 2) Draw machine effects: psi_j ~ N(1,1)
preserve
clear
set obs `J'
gen machine = _n
gen psi     = rnormal(1,1)
tempfile machines
save `machines'
restore


* 3) Draw period effects: gamma_t ~ N(0,0.5)
preserve
clear
set obs `T'
gen period = _n
gen gamma  = rnormal(0, sqrt(0.5))
tempfile periods
save `periods'
restore

* 4) Build worker-period panel

clear
set obs `=`N'*`T''
gen worker = ceil(_n/`T')
gen period = mod(_n-1, `T') + 1

merge m:1 worker using `workers', nogen
merge m:1 period using `periods', nogen


* 5) Random machine assignment:
gen machine = ceil(`J' * runiform())

merge m:1 machine using `machines', nogen


* 6) Generate noise 
gen eps = rnormal(0, sqrt(0.2))


* 7) Generate outcome

gen y = alpha + psi + gamma + eps


*** Estimate FE model ***

reg y i.worker i.machine i.period, vce(cluster worker)

* Recover worker effects from the regression
preserve
    keep worker alpha
    bys worker: keep if _n == 1
    tempfile true_workers
    save `true_workers'
restore

clear
set obs `N'
gen worker = _n
gen alpha_hat = 0

forvalues i = 2/`N' {
    replace alpha_hat = _b[`i'.worker] if worker == `i'
}

merge 1:1 worker using `true_workers', nogen

* Center both so the arbitrary normalization does not matter
summ alpha, meanonly
gen alpha_true_c = alpha - r(mean)

summ alpha_hat, meanonly
gen alpha_hat_c = alpha_hat - r(mean)

twoway ///
    (scatter alpha_hat_c alpha_true_c) ///
    (function y = x, range(alpha_true_c)), ///
    xtitle("True alpha_i (centered)") ///
    ytitle("Estimated alpha_i (centered)")
	
graph export "figures/scatter_1.pdf", replace

 * Recover machine effects from the regression
use `machines', clear

gen psi_hat = 0
forvalues j = 2/`J' {
    replace psi_hat = _b[`j'.machine] if machine == `j'
}

summ psi, meanonly
gen psi_true_c = psi - r(mean)

summ psi_hat, meanonly
gen psi_hat_c = psi_hat - r(mean)

twoway ///
    (scatter psi_hat_c psi_true_c) ///
	(function y = x, range(psi_true_c)), ///
    xtitle("True psi_j (centered)") ///
    ytitle("Estimated psi_j (centered)")
	
graph export "figures/scatter_2.pdf", replace

 

 *** 8 *** 
 
 
 
 **** Simulate data ****

clear all                 
set more off 

set seed 240505

* Dimensions
local N = 10
local T = 5
local J = 20

* 1) Draw worker effects: alpha_i ~ N(0,1)
preserve
clear
set obs `N'
gen worker = _n
gen alpha  = rnormal(0,1)
tempfile workers
save `workers'
restore

* 2) Draw machine effects: psi_j ~ N(1,1)
preserve
clear
set obs `J'
gen machine = _n
gen psi     = rnormal(1,1)
tempfile machines
save `machines'
restore


* 3) Draw period effects: gamma_t ~ N(0,0.5)
preserve
clear
set obs `T'
gen period = _n
gen gamma  = rnormal(0, sqrt(0.5))
tempfile periods
save `periods'
restore

* 4) Build worker-period panel

clear
set obs `=`N'*`T''
gen worker = ceil(_n/`T')
gen period = mod(_n-1, `T') + 1

merge m:1 worker using `workers', nogen
merge m:1 period using `periods', nogen


* 5) Random machine assignment:
gen machine = ceil(`J' * runiform())

merge m:1 machine using `machines', nogen


* 6) Generate noise 
gen eps = rnormal(0, sqrt(0.2))


* 7) Generate outcome

gen y = alpha + psi + gamma + eps


*** Estimate FE model ***

reg y i.worker i.machine i.period, vce(cluster worker)

* Recover worker effects from the regression
preserve
    keep worker alpha
    bys worker: keep if _n == 1
    tempfile true_workers
    save `true_workers'
restore

clear
set obs `N'
gen worker = _n
gen alpha_hat = 0

forvalues i = 2/`N' {
    replace alpha_hat = _b[`i'.worker] if worker == `i'
}

merge 1:1 worker using `true_workers', nogen

* Center both so the arbitrary normalization does not matter
summ alpha, meanonly
gen alpha_true_c = alpha - r(mean)

summ alpha_hat, meanonly
gen alpha_hat_c = alpha_hat - r(mean)

twoway ///
    (scatter alpha_hat_c alpha_true_c) ///
    (function y = x, range(alpha_true_c)), ///
    xtitle("True alpha_i (centered)") ///
    ytitle("Estimated alpha_i (centered)")
	
graph export "figures/scatter_3.pdf", replace

 * Recover machine effects from the regression
use `machines', clear

gen psi_hat = 0
forvalues j = 2/`J' {
    replace psi_hat = _b[`j'.machine] if machine == `j'
}

summ psi, meanonly
gen psi_true_c = psi - r(mean)

summ psi_hat, meanonly
gen psi_hat_c = psi_hat - r(mean)

twoway ///
    (scatter psi_hat_c psi_true_c) ///
	(function y = x, range(psi_true_c)), ///
    xtitle("True psi_j (centered)") ///
    ytitle("Estimated psi_j (centered)")
 
graph export "figures/scatter_4.pdf", replace
 









