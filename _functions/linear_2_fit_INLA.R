################################################################################
# This script contains functions for data management and fitting linear model
# with INLA package. The functions include:
#
# inclusion_IgGmonth;
# sort_colonization;
# prep_dataINLA;
# fit_INLA; and
# pluck_estINLA.
################################################################################

# Function for including subjects based on exposure measurement:
# whether IgG measurement at the selected month (7m or 13m) is available
inclusion_IgGmonth <- function(Isratrid_Export, IgG_month){
  
  # Check data on IgG measurement and log-transform IgG concentration
  df_aby <- Isratrid_Export %>% 
    select(-c(Pnc_POS_NEG, First_serotype_detected, Second_serotype_detected, Third_serotype_detected)) %>% 
    pivot_longer(Concentration_serotype_1:Concentration_serotype_23F, 
                 names_to = "Serotype", values_to = "Concentration") %>% 
    separate(Serotype, into = c("c", "s", "Serotype"), sep = "_") %>% 
    select(-c(c,s)) %>% 
    mutate(Serotype = as.factor(Serotype),
           log_GMC = log(Concentration))
  
  # Focus on subjects with IgG measurement at the selected month (7m or 13m)
  df_aby <- df_aby %>% filter(Age_in_months == IgG_month)
  vec_id <- unique(df_aby$Study_key)
  
  # Get colonization for included subjects
  df_colonization <- Isratrid_Export %>%
    filter(Study_key %in% vec_id) %>% # Checked: ok to filter here, doesn't change included case base
    select(-c(Visit_number:Visit_Month, Third_serotype_detected, Concentration_serotype_1:Concentration_serotype_23F))
  
  ls_res <- list(df_aby, df_colonization)
  names(ls_res) <- c("df_aby", "df_colonization")
  return(ls_res)
}
################################################################################
# End of function
################################################################################


# Function for getting colonization observations and censoring subsequent same-
# serotype colonization such that only new acquisitions are included 
sort_colonization <- function(df_colonization, subject_id, vec_st, df_aby){
  
  # Create negation function
  `%notin%` <- Negate(`%in%`)
  
  # Focus on one study subject, on PCV7 serotypes
  x <- df_colonization %>% filter(Study_key == subject_id) 
  vec_age <- unique(x$Age_in_months)
  x <- x %>% 
    pivot_wider(names_from = "Age_in_months", values_from = "Pnc_POS_NEG") %>% 
    pivot_longer(First_serotype_detected:Second_serotype_detected, names_to = "Detection", values_to = "Serotype")
  x1 <- x %>% filter(Serotype %in% vec_st) %>% select(-Detection)
  
  # Determine total no. of follow-up visits attended
  n_col <- ncol(x1)
  n_visit <- n_col - 5
  
  if(nrow(x1) == 0){
    vec_st_neg <- vec_st
    x2 <- NULL
  } else{
    vec_st_neg <- vec_st[vec_st %notin% unique(x1$Serotype)]
    # Recode POS & NEG
    x2 <- x1 %>% select(1:4, Serotype, 5:n_col-1) %>% 
      pivot_longer(6:n_col, names_to = "Age_in_months", values_to = "Colonization") %>% 
      mutate(Colonization = ifelse(Colonization=="POS", 1, 0),
             Age_in_months = as.numeric(Age_in_months))
    x2$Colonization <- replace(x2$Colonization, is.na(x2$Colonization), 0)
  }
  
  # Fill in negative results for undetected PCV7 serotypes
  id <- x[1,1:4]
  y1 <- rbind(id, id[rep(1, n_visit-1), ])
  y2 <- data.frame(Age_in_months = vec_age)
  y3 <- data.frame(Colonization = rep(0, n_visit))
  
  ls_neg <- NULL
  for(i in 1: length(vec_st_neg)){
    z <- data.frame(Serotype = rep(vec_st_neg[i], n_visit))
    ls_neg[[i]] <- cbind(y1, z, y2, y3) 
  }
  df_res <- bind_rows(ls_neg, x2)
  
  # Once colonization by a serotype occurs, remove subsequent observations for that serotype 
  df_res <- df_res %>%
    group_by(Serotype) %>% 
    filter(row_number() <= which.max(Colonization)|all(Colonization==0))
  
  # Combine with IgG measurement
  df_IgG <- df_aby %>%
    filter(Study_key == subject_id, Serotype %in% vec_st) %>%
    select(Study_key, Serotype, log_GMC)
  df_res <- df_res %>% left_join(df_IgG, by = c("Study_key", "Serotype"))
  
  # Return output data frame
  return(df_res)
}
################################################################################
# End of function
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
