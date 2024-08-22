################################################################################
# Function for extracting estimates from INLA model
################################################################################
pluck_estINLA <- function(modINLA, M = c(1:4), vec_serotype){
  # Set model description
  if(M == 1) { mod <- "Ethnic,logGMC" }
  if(M == 2) { mod <- "Ethnic,logGMC,agec" }
  if(M == 3) { mod <- "Ethnic,logGMC,st_ethnic" }
  if(M == 4) { mod <- "Ethnic,logGMC,st_ethnic,agec" }
  
  if(M == 1 | M == 2){
    mean_logGMC <- modINLA[["summary.fixed"]]["log_GMC", "mean"]
    sd_logGMC <- modINLA[["summary.fixed"]]["log_GMC", "sd"]
    
    res <- modINLA$summary.random$Serotype2 %>%
      mutate(mod = mod,
             mean_b = mean + mean_logGMC,
             var = sd^2 + sd_logGMC^2,
             upper = mean_b + 1.96*sqrt(var),
             lower = mean_b - 1.96*sqrt(var)) %>%
      cbind(vec_serotype) %>% 
      select(mod, vec_serotype, mean_b, var, upper, lower)
    names(res)[2] <- c("serotype")
  }
  
  if(M == 3 | M == 4){
    mean_logGMC <- modINLA[["summary.fixed"]]["log_GMC", "mean"]
    sd_logGMC <- modINLA[["summary.fixed"]]["log_GMC", "sd"]
    
    beta_serotype <- as.data.frame(modINLA$summary.random$Serotype2) %>% 
      select(ID, mean, sd)
    beta_serotype_ethnicity <- as.data.frame(modINLA$summary.random$Serotype_Ethnicity2) %>% 
      select(ID, mean, sd) %>% 
      mutate(ID2 = ID) %>% 
      separate(ID, into = c("serotype", "ethnicity"), sep = "_")
    
    res <- left_join(beta_serotype_ethnicity, beta_serotype, by = c("serotype" = "ID")) %>% 
      mutate(mod = mod,
             mean_b = mean.x + mean.y + mean_logGMC,
             var = sd.x^2 + sd.y^2 + sd_logGMC^2,
             upper = mean_b + 1.96*sqrt(var),
             lower = mean_b - 1.96*sqrt(var)) %>% 
      select(mod, serotype, ethnicity, mean_b, var, upper, lower)
  }
  
  return(res)
}
################################################################################
# End of function
################################################################################
