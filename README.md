# Estimating serotype-specific association between the concentration of serum antibodies and risk of pneumococcal colonization
This project aims to estimate the serotype-specific relative risk (RR) of colonization comparing a higher-valency pneumococcal conjugate vaccine (PCV2) to a lower-valency one (PCV1). The relative vaccine effectiveness (VE) against colonization can be obtained by 1-RR. 

For example, assuming that PCV7 confers 60% (RR=0.4) protection against colonization by a serotype, an estimated RR of 1.10 when comparing PCV2 with PCV1 can be interpreted as 10% reduction in the absolute VE for PCV2, i.e., PCV2 has an absolute VE of 56% (1-0.4*1.1).


We formulated two hierarchical Bayesian modelling approaches:

(1) Linear model

log(RR) = b1 * log(GMR)

where 
* b1 is the serotype-specific, ethnicity-specific coefficient; and 
* log(GMR) is the natural-logarithm of the Geometric Mean Concentration Ratio between PCV2 and PCV1. The GMCs give the antibody levels following immunization with PCV2 and PCV1 reported in a head-to-head clinical trial.

(2) Change point model

Extending from the linear model, we assumed that for a given serotype, the risk of colonization remained constant until serum IgG reached a threshold concentration (the change point). After the change point, the risk of colonization decreased with increasing concentration of serum IgG, following a log-log linear relationship.


## Data

* To fit the hierarchical Bayesian model, this project uses the data from [Dagan et al. 2013](https://pubmed.ncbi.nlm.nih.gov/23804191/), [Dagan et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27422342/).
* To estimate the RR for PCV13 vs. PCV7, PCV15 vs. PCV13, and PCV20 vs. PCV13, this study uses the summary-level immunogenicity data from head-to-head clinical trials, identified from a systematic review [Fitch et al. 2024](https://elischolar.library.yale.edu/ysphtdl/2390/). The data were curated from https://wisspar.com, and https://classic.clinicaltrials.gov.


## Steps

Here are the steps for the two approaches:

(1) Linear model
* Getting logGMR from trial data
  * 1_get_wtGMR.Rmd
* Fitting hierarchical Bayesian model (approximation with INLA) to obtain b1 
  * 2_fit_INLA.Rmd
* Calculating RR
  * 3_get_rVE.Rmd
  
(2) Change point model
* Fitting hierarchical Bayesian model (using rjags) to obtain b1 and change point
  * 1_fit_rjags.Rmd
* Calculating RR
  * 2_get_rVE.Rmd
