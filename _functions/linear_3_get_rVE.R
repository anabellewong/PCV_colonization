################################################################################
# Function for calculating relative VE from poled GMR
################################################################################
get_rVE <- function(est, df_GMR, vec_st, vax_dose = c("postprim", "postboost")){
  df_rVE <- left_join(est, df_GMR, by = c("serotype" = "Serotype")) %>% 
    mutate(serotype = fct_relevel(as.factor(serotype), vec_st),
           data = fct_relevel(as.factor(data), c("post-primary, obs 7-24m",
                                                 #"post-primary, obs 13-24m",
                                                 "post-booster, obs 13-24m")),
           LogRR = mean_b * LogGMR,
           se_LogRR = sqrt(var*var_LogGMR + var_LogGMR*mean_b^2 + var*LogGMR^2),
           uci_LogRR = LogRR + 1.96*se_LogRR,
           lci_LogRR = LogRR - 1.96*se_LogRR)
  if(vax_dose == "postprim"){
    df_rVE <- df_rVE %>% filter(data != "post-booster, obs 13-24m")
  }
  else{
    df_rVE <- df_rVE %>% filter(data == "post-booster, obs 13-24m")
  }
  return(df_rVE)
}
################################################################################
# End of function
################################################################################
