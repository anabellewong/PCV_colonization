################################################################################
# Function for preparing data for fitting INLA models and cutting age into age 
# groups; the cuts depend on whether exposure if IgG measured at 7m or 13m
################################################################################
prep_dataINLA <- function(df_colonization, IgG_month){
  df_INLA <- df_colonization %>%
    mutate(log_GMC2=log_GMC, # just a copy since INLA doesn't allow doing the same var twice
           log_GMC3=log_GMC, # just a copy since INLA doesn't allow doing the same var twice
           Serotype=as.factor(Serotype),
           Serotype2 = Serotype, # just a copy
           Ethnicity=as.factor(Ethnicity),
           Serotype_Ethnicity= as.factor(paste(Serotype, Ethnicity,sep='_')),
           Serotype_Ethnicity2 = Serotype_Ethnicity) # just a copy

    if(IgG_month == 7){
      df_INLA <- df_INLA %>%
        mutate(agec = case_when(Age_in_months>=7 & Age_in_months<12 ~ 1,
                                Age_in_months>=12 & Age_in_months<13 ~ 2,
                                Age_in_months>=13 & Age_in_months<18 ~ 3,
                                Age_in_months>=18 & Age_in_months<Inf ~ 4),
               agec = as.factor(agec))
    }
    else{
      df_INLA <- df_INLA %>%
        mutate(agec = case_when(Age_in_months>=13 & Age_in_months<16 ~ 1,
                                Age_in_months>=16 & Age_in_months<19 ~ 2,
                                Age_in_months>=19 & Age_in_months<22 ~ 3,
                                Age_in_months>=22 ~ 4),
               
               agec = as.factor(agec))
    }

  stopifnot(sum(is.na(df_INLA$agec))==0)
  return(df_INLA)
}
################################################################################
# End of function
################################################################################