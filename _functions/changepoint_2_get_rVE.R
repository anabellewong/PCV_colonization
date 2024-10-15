################################################################################
# This script contains functions for getting relative VE by pooling RR across 
# studies with inverse-variance weights. The RR came from comparing the risk of 
# colonization for different PCVs predicted by the change point model. The 
# functions include:
#
# log_transform_igg;
# predict_col;
# calc_RR; and
# get_wtRR.
################################################################################

# Log-transform IgG levels
log_transform_igg <- function(df_postboost, df_map){
  df_postboost$lower_limit[which(df_postboost$lower_limit==0)] <- 0.01
  df_postboost$upper_limit[which(df_postboost$upper_limit==0)] <- 0.01
  df_postboost$GMC[which(df_postboost$GMC==0)] <- 0.01
  df_postboost <- df_postboost %>% 
    mutate(across(c(GMC:lower_limit), function(x) log(x)),
           se_logGMC = abs(upper_limit - lower_limit)/(1.96*2)) %>% 
    left_join(df_map, by = "serotype")
  return(df_postboost)
}
################################################################################
# End of function
################################################################################


# Get predicted risks from rjags model, which contains posteriors of
# b0 (13) starting with index: 1
# b1 (13) starting with index: n_serotype + 1
# b3 (13) starting with index: 2*n_serotype + 1
# cp (13) starting with index: 3*n_serotype + 1

predict_col <- function(mod, st.index, mean_logGMC, sd_logGMC, ethnic = c("jewish", "bedouin")){
  # Take one serotype, get the 30000 posterior for b0
  vec_b0 <- mod$posterior_samples.all[ , st.index]
  vec_logGMC <- rnorm(length(vec_b0), mean = mean_logGMC, sd = sd_logGMC)
  
  # Get b3
  vec_b3 <- mod$posterior_samples.all[ ,2*n_serotype + st.index]
  
  # Get b2 (b2 = b1*logGMC for GMC after change point)
  vec_b1 <- mod$posterior_samples.all[ ,n_serotype + st.index]
  delta_logGMC <- vec_logGMC - mod$posterior_samples.all[ , 3*n_serotype + st.index]
  vec_b2 <- vec_b1 * delta_logGMC
  
  # Is the GMC above change point?
  vec_cp <- vec_logGMC > mod$posterior_samples.all[ , 3*n_serotype + st.index] 
  
  # Get the log_pred_col
  if(ethnic == "bedouin"){
    log_pred_col <- vec_b2 * vec_cp + (vec_b0+vec_b3)
  }
  if(ethnic == "jewish"){
    log_pred_col <- vec_b2 * vec_cp + (vec_b0)
  }
  
  return(log_pred_col)
}
################################################################################
# End of function
################################################################################


# Calculate RR relying on the function predict_col() and the rjags model object
calc_RR <- function(df_postboost, mod, ethnic = c("jewish", "bedouin")){
  
  # Create an empty list
  ls_RR <- setNames(vector(mode = "list", length = n_serotype), nm = c(1:13))
  
  set.seed(102405L)
  for(i in 1:n_serotype){
    data <- df_postboost %>% filter(st.index == i)
    vec_study <- unique(data$study_id)
    
    for(nm in vec_study){
      d <- data %>% filter(study_id == nm)
      d1 <- d %>% filter(vaccine == PCV1) # lower-valency
      d2 <- d %>% filter(vaccine == PCV2) # higher-valency
      
      # Predict log-risk from model object, choosing the ethnic-specific intercept
      if(ethnic == "bedouin"){
        log_risk_1 <- predict_col(mod, st.index = i, mean_logGMC = d1$GMC, sd_logGMC = d1$se_logGMC, ethnic = "bedouin")
        log_risk_2 <- predict_col(mod, st.index = i, mean_logGMC = d2$GMC, sd_logGMC = d2$se_logGMC, ethnic = "bedouin")
      }
      if(ethnic == "jewish"){
        log_risk_1 <- predict_col(mod, st.index = i, mean_logGMC = d1$GMC, sd_logGMC = d1$se_logGMC, ethnic = "jewish")
        log_risk_2 <- predict_col(mod, st.index = i, mean_logGMC = d2$GMC, sd_logGMC = d2$se_logGMC, ethnic = "jewish")
      }
      
      # Calculate RR comparing higher-valency vs. lower-valency PCV
      RR <- log_risk_2 - log_risk_1
      
      # Get mean & 95ci from 30000 iter
      vec_res <- c(quantile(RR, c(0.5, 0.025, 0.975)), sd(RR))
      names(vec_res) <- c("mean", "lci", "uci", "sd")
      ls_RR[[i]][[nm]] <- vec_res
    }
    ls_RR[[i]] <- bind_rows(ls_RR[[i]], .id = "study_id")
  }
  df_RR <- bind_rows(ls_RR, .id = "st.index") %>% 
    mutate(st.index = as.numeric(st.index))
  return(df_RR)
}
################################################################################
# End of function
################################################################################


# Function for calculating inverse-variance weighted average across studies
get_wtRR <- function(df_RR, Comparison){
  df_wtRR <- df_RR %>% 
    mutate(st.index = as.numeric(st.index),
           inv_var = 1/(sd^2)) %>% 
    group_by(st.index) %>% 
    mutate(total_inv_var = sum(inv_var),
           wt_inv = inv_var/ total_inv_var,
           wt_RR = mean*wt_inv) %>% 
    summarise(sum_wt_RR = sum(wt_RR),
              var_wt_RR = 1/total_inv_var,
              se_wt_RR = sqrt(var_wt_RR),
              wt_RR_uci = sum_wt_RR + 1.96*se_wt_RR,
              wt_RR_lci = sum_wt_RR - 1.96*se_wt_RR) %>% 
    unique() %>% 
    left_join(df_map, by = "st.index") %>% 
    mutate(serotype = as.factor(serotype),
           serotype = fct_relevel(serotype, vec_serotype)) %>% 
    select(serotype, sum_wt_RR, wt_RR_lci, wt_RR_uci, se_wt_RR)
  df_wtRR <- df_wtRR[,-1]
  names(df_wtRR) <- c("serotype", "LogRR", "lci_LogRR", "uci_LogRR", "se_LogRR")
  df_wtRR <- df_wtRR %>% mutate(Comparison = Comparison)
  return(df_wtRR)
}
################################################################################
# End of function
################################################################################
