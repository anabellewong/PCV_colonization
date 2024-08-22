################################################################################
# Function for fitting INLA model
################################################################################
fit_INLA <- function(dataINLA, mod = c(1:4)){
  if(mod == 1){
    f <- as.formula("Colonization ~ Ethnicity + log_GMC +
                  f(Serotype, model = 'iid', hyper = prior_hyper) +
                  f(Serotype2, log_GMC, model = 'iid', hyper = prior_hyper)")
  }
  if(mod == 2){
    f <- as.formula("Colonization ~ Ethnicity + log_GMC + agec + 
                      f(Serotype, model = 'iid', hyper = prior_hyper) +
                      f(Serotype2, log_GMC, model = 'iid', hyper = prior_hyper)")
  }
  if(mod == 3){
    f <- as.formula("Colonization ~ Ethnicity + log_GMC +
                      f(Serotype, model = 'iid', hyper = prior_hyper) +
                      f(Serotype2, log_GMC, model = 'iid', hyper = prior_hyper) +
                      f(Serotype_Ethnicity, model = 'iid', hyper = prior_hyper) +
                      f(Serotype_Ethnicity2, log_GMC, model = 'iid', hyper = prior_hyper)")
  }
  if(mod == 4){
    f <- as.formula("Colonization ~ Ethnicity + log_GMC + agec + 
                      f(Serotype, model = 'iid', hyper = prior_hyper) +
                      f(Serotype2, log_GMC, model = 'iid', hyper = prior_hyper) +
                      f(Serotype_Ethnicity, model = 'iid', hyper = prior_hyper) +
                      f(Serotype_Ethnicity2, log_GMC, model = 'iid', hyper = prior_hyper)")
  }
  mod <- inla(formula = f, 
                  data = dataINLA,
                  family="binomial", 
                  control.family=list(link='log'),
                  control.fixed = prior_fixed,
                  control.compute = list(waic = TRUE))
  return(mod)
}
################################################################################
# End of function
################################################################################
