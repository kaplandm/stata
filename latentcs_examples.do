* Author: David M. Kaplan
* Updated: 17aug2022
* Purpose: examples of latentcs.ado
* 1. Example from  help latentcs
* 2. Newborn health (one-minute apgar)
* 3. Simulated data
* 4. Mental health


* To install:
*net from https://kaplandm.github.io/stata
*net describe latentcs
*net install latentcs
*net get latentcs


* Specify directory w/ .ado (if not already in a standard directory; use command
*    sysdir   to check standard directories.
*cd "C:\Users\kaplandm\code"

*.tex output:  sjlog do latentcs_examples , clear replace

// local SAVEFILES 1 // 1=save graphs; 0=do not
set more off
set linesize 67
set level 95


*********
* 1. Example from   help latentcs
*********
set scheme sj
webuse hanley , clear
describe disease rating
latentcs rating , by(disease)

*********
* 2. Example: newborn health (one-minute apgar score)
*********
use "http://fmwww.bc.edu/ec-p/data/wooldridge/bwght2", clear
gen smoke = (cigs>0)
latentcs omaps , by(smoke) switch level(90)
latentcs omaps , by(smoke) switch level(80)
gen omapsagg = omaps - 2
replace omapsagg = 1 if omapsagg<=1
replace omapsagg = 5 if omapsagg>=5
latentcs omapsagg , by(smoke) switch level(90)
capture drop omapsagg2
gen omapsagg2 = (omaps>5)
latentcs omapsagg2 , by(smoke) switch level(90)

*********
* 3. Example with simulated data
*********
* Set RNG seed for replication
set seed 112358

* Set up DGP, dataset
clear
local nt 100
local nc 100
disp `=`nc'+`nt''
set obs `=`nc'+`nt''
gen treated = 0
replace treated=1 if _n<=`nt'

* Generate outcomes
* Case 1: no change in distribution
gen N01 = rnormal(0,1)
gen y1 = 1+(N01>-1.5)+(N01>-1)+(N01>-0.5)+(N01>0)+(N01>0.5)+(N01>1)+(N01>1.5)
* Case 2: treatment effect above median
replace N01 = N01+2 if treated==1 & N01>0
gen y2 = 1+(N01>-1.5)+(N01>-1)+(N01>-0.5)+(N01>0)+(N01>0.5)+(N01>1)+(N01>1.5)

* Run latentcs
latentcs y1 , by(treated) // no actual effect (but alpha probability of error)
latentcs y2 , by(treated) // effect only above median (0)


*********
* 4. Example from paper: mental health
*********
set more off
import delimited "https://docs.google.com/uc?id=11bfNXYAhuhye9ICA0O6BL7voHgOEOvht&export=download" , varnames(1) clear
gen female = (sex==2)
drop sex
gen black = (racenew==20)
gen educ2 = educ>302
replace educ2=. if educ==0 | educ>=996
gen K6 = 25-(aeffort+ahopeless+anervous+arestless+asad+aworthless)
replace K6=. if aeffort>4 | ahopeless>4 | anervous>4 | arestless>4 | asad>4 | aworthless>4
gen poverty = (pooryn==2)
replace poverty=. if pooryn==9
replace health=. if health>5
replace health = 6-health
* Now: female=TRUE if female
*      black=1 if racenew=20=Black (o/w 0 or .)
*      educ2=1 if any postsecondary (o/w 0 or .)
*      poverty=1 if below poverty threshold (o/w 0 or .)
*      K6 = reversed Kessler-6 for nonspecific psychological distress
*           from 1=worst to 25=best (o/w .)
*           1-12: maybe serious mental illness
*           13-20: moderate mental distress
*           21-25: probably no mental health issues
*      health = {1=poor,2=fair,3=good,4=very good,5=excellent} (o/w .)
* (and year, nhispid, sampweight [for K6], perweight [for health] as before)

preserve
drop if missing(K6) | missing(poverty) | female
tab year poverty if !missing(K6) , missing
set level 90
* Weighted
latentcs K6 [aw=sampweight] if  poverty & !female , by(year) switch
latentcs K6 [aw=sampweight] if !poverty & !female , by(year) switch
tab K6 year if poverty & !female [aw=sampweight], missing column nofreq
* Single category (18)
gen K6gt18 = (K6>18)
latentcs K6gt18 [aw=sampweight] if  poverty & !female , by(year) switch ndraw(100000) noplot
* Aggregate small/similar categories
gen K6agg1 = 1+floor((K6-1)/5)
latentcs K6agg1 [aw=sampweight] if  poverty & !female , by(year) switch noplot
* Different aggregation
gen K6agg2 = K6-11
replace K6agg2=1 if K6agg2<1
replace K6agg2=8 if K6agg2>8 & !missing(K6)
latentcs K6agg2 [aw=sampweight] if  poverty & !female , by(year) switch ndraw(100000) noplot
restore

* End of file
