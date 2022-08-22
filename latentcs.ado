*! version 0.1 12aug2022
program latentcs, rclass sortpreserve byable(recall)
  version 11
  syntax varname(numeric default=none) [if] [in] [aweight fweight] , BY(varname) [SWitch Level(cilevel) SEED(integer 112358) NDRAW(integer 10000) noPLOT ]
  * Basic setup and check for problems before calling mata function
  if ("`level'"=="") {
    local level = c(level)
    if ("`level'"=="") {
      local level = 90
    }
  }
  if (`level'<10 | `level'>99.99) { // matches regress
    di as error "{bf:level(`level')} is out of range: must be between 10 and 99.99 inclusive"
    exit 198
  }
  marksample touse
  qui replace `touse'=0 if missing(`by')
  qui replace `touse'=0 if missing(`varlist')
  qui count if `touse'
  if !r(N) { //taken from ksmirnov.ado
    error 2000
  }
  qui tab `by' if `touse'
  if r(r) != 2 {
    di as error "`by' takes on " r(r) " values, not 2"
    exit 450
  }

  * Check if `by' is string or numeric; tab , matrow() not supported for string
  local grpvar "`by'"
  capture confirm string var `by'
  if (_rc==0) { //string
    tempname bynum
    encode `by' if `touse' , g(`bynum')
    local grpvar "`bynum'"
  }

  * Compute ECDFs, etc.
  tempname mfreq n0n1 epmfs ecdfs ycats vals
  qui tab `varlist' `grpvar' if `touse' [`weight'`exp'] , matcell(`mfreq') matrow(`ycats') matcol(`vals')
  local v0 = `vals'[1,1]
  local v1 = `vals'[1,2]
  if ("`switch'"!="") {
    matrix `mfreq' = `mfreq'[1..rowsof(`mfreq'),2] , `mfreq'[1..rowsof(`mfreq'),1]
    local v0 = `vals'[1,2]
    local v1 = `vals'[1,1]
  }
  matrix `n0n1' = J(1,rowsof(`mfreq'),1)*`mfreq'
  matrix `epmfs' = `mfreq'[1..rowsof(`mfreq'),1]/`n0n1'[1,1] , `mfreq'[1..rowsof(`mfreq'),2]/`n0n1'[1,2]
  * Unless fweight, use subsample counts for subsample sizes (instead of weighted counts)
  if ("`weight'"!="fweight") {
    count if `grpvar'==`v0' & `touse'
    matrix `n0n1'[1,1] = r(N)
    count if `grpvar'==`v1' & `touse'
    matrix `n0n1'[1,2] = r(N)
  }
  local grp0lab : label `grpvar' `v0'
  local grp1lab : label `grpvar' `v1'
  * Move estimated PMFs to Mata and compute corresponding CDFs
  mata: epmfs = st_matrix("`epmfs'")
  mata: ecdfs = runningsum(epmfs[.,1]) , runningsum(epmfs[.,2])
  mata: st_matrix( "`ecdfs'" , ecdfs )

  * Main subroutine call
  tempname seedname retest retcs
  local `seedname' = c(seed)
  mata: latentcsmain("`grpvar'", `level', "`touse'", "`n0n1'", epmfs, ecdfs, `seed', `ndraw', "`retest'", "`retcs'")
  set seed ``seedname''

  * Store r() values
  return scalar level = `level'
  tempname nmat
  matrix `nmat' = (`n0n1'[1,1]+`n0n1'[1,2] , `n0n1'[1,1] , `n0n1'[1,2])
  matrix colnames `nmat' = Ntotal Ngroup1 Ngroup2
  return matrix N = `nmat'

  * Display results
  matrix colnames `retest' = from to
  matrix colnames `retcs' = from to
  di as text "Results for quantile ranges where distribution for"
  di as text "  `by'=`grp1lab' has higher tau-quantile of `varlist' than `by'=`grp0lab'"
  di as text " "
  if (`retest'[1,1]==.) {
    di as text "Estimated: none"
  }
  else {
    di as text "Estimated:"
    matlist `retest' , colorcoleq(res) lines(none) names(columns) noblank
  }
  di as text " "
  if (`retcs'[1,1]==.) {
    di as text "Inner `level'% confidence set: empty"
  }
  else {
    di as text "Inner `level'% confidence set:"
    matlist `retcs' , colorcoleq(res) lines(none) names(columns) noblank
  }

  * Save in r()
  return matrix tau_bw_est = `retest'
  return matrix tau_bw_cs = `retcs'

  * Plot (or not)
  if ("`plot'"!="noplot") {
    tempname ycatvar ecdf0var ecdf1var
    qui gen `ycatvar'  = `ycats'[_n,1] in 1/`=rowsof(`ycats')'
    qui gen `ecdf0var' = `ecdfs'[_n,1] in 1/`=rowsof(`ycats')'
    qui gen `ecdf1var' = `ecdfs'[_n,2] in 1/`=rowsof(`ycats')'
    graph twoway (line `ecdf0var' `ycatvar' in 1/`=rowsof(`ycats')' , connect(stairstep) ) ///
                 (line `ecdf1var' `ycatvar' in 1/`=rowsof(`ycats')' , connect(stairstep) ) ///
          , legend(lab(1 "`by'=`grp0lab'") lab(2 "`by'=`grp1lab'")) xtitle("`varlist'") ytitle("CDF")
  }
end
*
*
*
mata:
void latentcsmain(string scalar grpname, real scalar conflevel, string scalar touse, string scalar n0n1name, real matrix epmfs, real matrix ecdfs, real scalar seed, real scalar ndraw, string scalar retestname, string scalar retcsname) {
  real scalar i, j, k, eps, JJ, n0, n1, alpha, talpha, tbeta // counters, sample sizes, etc.
  real matrix tmp, n0n1, sigma0hat, sigma1hat, cs, est, ecdfsraw
  // real matrix epmfs, ecdfs
  string colvector grp // from dataset; only used if string
  real colvector F0hat, F1hat, conf0, conf1, ind // ECDFs, confidence limits, etc.

  // Load from Stata
  n0n1 = st_matrix(n0n1name)
  n0 = n0n1[1]
  n1 = n0n1[2]
  alpha = (100 - conflevel) / 100

  // Aggregate estimated intervals
  est = J(0,2,.)
  JJ = rows(ecdfs)
  F0hat = ecdfs[.,1]
  F1hat = ecdfs[.,2]
  for (j=1; j<=(JJ-1); j++) {
    if (F1hat[j]<F0hat[j]) {
      if (rows(est)>0) {
        if (F1hat[j]<est[rows(est),2]) {
          est[rows(est),2] = F0hat[j]
        }
        else {
          est = est \ (F1hat[j], F0hat[j])
        }
      }
      else {
        est = est \ (F1hat[j], F0hat[j])
      }
    }
  }
  if (rows(est)>0) {
    st_matrix(retestname, est)
  }
  else {
    st_matrix(retestname, J(1,2,.))
  }

  // Remove categories with zero observations (in one group or another)
  ecdfsraw = ecdfs
  ecdfs = J(0,2,.)
  for (i=1; i<=rows(ecdfsraw); i++) {
    if (epmfs[i,1]!=0) { //otherwise, just skip
      if (epmfs[i,2]!=0 | i==1) {
        ecdfs = ecdfs \ ecdfsraw[i,.]
      }
      else {
        if (epmfs[i-1,1]==0) {
          ecdfs = ecdfs \ ecdfsraw[i,.]
        }
        else {
          ecdfs[rows(ecdfs),.] = ecdfsraw[i,.]
        }
      }
    }
  }
  if (ecdfs[rows(ecdfs),1]<(1-2*epsilon(1)) | ecdfs[rows(ecdfs),2]<(1-2*epsilon(1))) {
    ecdfs = ecdfs \ (1,1)
  }

  // Finite-sample adjustment: if applicable, change 0 to 1/n, 1 to 1-1/n
  k = max(ecdfs[1,.])
  ecdfs = ecdfs + ( (ecdfs[.,1]:==0)*min((k,1/n0)) , (ecdfs[.,2]:<=(2*epsilon(1)))*min((k,1/n1)) )
  k = min(ecdfs[rows(ecdfs)-1,.])
  if (k==1) {
    k = 1 - 1/max((n0,n1))
  }
  ecdfs = ecdfs - ( (ecdfs[.,1]:>=(1-2*epsilon(1)))*min((1-k,1/n0)) , (ecdfs[.,2]:>=(1-2*epsilon(1)))*min((1-k,1/n1)) )
  JJ = rows(ecdfs)
  F0hat = ecdfs[.,1]
  F1hat = ecdfs[.,2]

  // Compute asymptotic covariance matrices
  sigma0hat = J(JJ-1, JJ-1, .)
  sigma1hat = sigma0hat
  for (i=1; i<=JJ-1; i++) {
    for (k=i; k<=JJ-1; k++) {
      sigma0hat[i,k] = F0hat[i]*(1-F0hat[k])
      sigma0hat[k,i] = sigma0hat[i,k]
      sigma1hat[i,k] = F1hat[i]*(1-F1hat[k])
      sigma1hat[k,i] = sigma1hat[i,k]
    }
  }

  // Correlation matrix (t-statistic covariance matrix)
  if (JJ>2) {
    sigma0t = corr(sigma0hat)
    sigma1t = corr(sigma1hat)
  }
  else {
    sigma0t = 1
    sigma1t = 1
  }

  // This should not happen anymore
  if (JJ>2 & (missing(cholesky(sigma0t)[1,1]) | missing(cholesky(sigma1t)[1,1]))) {
    stata(`"display as error "Covariance matrix not positive definite; can't compute inner CS""')
    printf("(Try aggregating categories to avoid categories with zero observations)\n")
    return
  }

  // Critical value simulation
  rseed(seed) // to be reproducible

  // Simulate/draw ndraw normal vectors (for first group t-statistics)
  if (JJ>2) {
    mvn0 = invnormal(uniform(ndraw,cols(sigma0t)))*cholesky(sigma0t)'  //per William Gould, p. 122 of https://doi.org/10.1177/1536867X0600600108
    PhiTmax = sort(normal(rowmax(mvn0)), 1)
    mvn0 = NULL
    // Compute \tilde\beta from paper (1-\tilde\beta = pointwise coverage prob)
    // Compute sample quantiles as in (1) in https://doi.org/10.1016/j.jeconom.2016.09.015
    k = floor(sqrt(1-alpha)*(ndraw+1))
    eps = (sqrt(1-alpha)*(ndraw+1))-k
    tbeta = 1 - ( (1-eps)*PhiTmax[k] + eps*PhiTmax[k+1] )
  }
  else {
    tbeta = 1 - sqrt(1-alpha)
  }
  // Confidence limits from Method 1, eqn (15)
  conf0 = J(JJ-1,1,.)
  for (i=1; i<=(JJ-1); i++) {
    conf0[i] = F0hat[i] - invnormal(1-tbeta)*sqrt(sigma0hat[i,i]/n0)
  }

  // Simulate/draw ndraw normal vectors (for second group t-statistics)
  if (JJ>2) {
    mvn1 = invnormal(uniform(ndraw,cols(sigma1t)))*cholesky(sigma1t)'  //per William Gould, p. 122 of https://doi.org/10.1177/1536867X0600600108
    PhiTmin = sort(normal(rowmin(mvn1)), 1)
    mvn1 = NULL
    // Compute \tilde\alpha from paper (1-\tilde\alpha = pointwise coverage prob)
    // Compute sample quantiles as in (1) in https://doi.org/10.1016/j.jeconom.2016.09.015
    k = floor((1-sqrt(1-alpha))*(ndraw+1))
    eps = ((1-sqrt(1-alpha))*(ndraw+1))-k
    talpha = (1-eps)*PhiTmin[k] + eps*PhiTmin[k+1]
  }
  else {
    talpha = 1-sqrt(1-alpha)
  }
  // Confidence limits from Method 1, eqn (15)
  conf1 = J(JJ-1,1,.)
  for (i=1; i<=(JJ-1); i++) {
    conf1[i] = F1hat[i] + invnormal(1-talpha)*sqrt(sigma1hat[i,i]/n1)
  }

  // Aggregate quantile index intervals into overall inner confidence set
  cs = J(0,2,.)
  for (j=1; j<=(JJ-1); j++) {
    if (conf1[j]<conf0[j]) {
      if (rows(cs)>0) {
        if (conf1[j]<cs[rows(cs),2]) {
          cs[rows(cs),2] = conf0[j]
        }
        else {
          cs = cs \ (conf1[j], conf0[j])
        }
      }
      else {
        cs = cs \ (conf1[j], conf0[j])
      }
    }
  }
  if (rows(cs)>0) {
    st_matrix(retcsname, cs)
  }
  else {
    st_matrix(retcsname, J(1,2,.))
  }
}
end // end mata
*EOF