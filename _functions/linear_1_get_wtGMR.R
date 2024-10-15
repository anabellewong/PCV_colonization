################################################################################
# This script contains functions for getting pooled estimates for Geometric 
# Mean Concentration Ratio (GMR). The functions include:
#
# transform_GMC;
# get_GMR;
# create_wt;
# get_wtGMR;
# patch_NA; and
# bind_df.
################################################################################

# Log-transform GMC and get se
transform_GMC <- function(GMC_data, vec_serotype){
  df <- GMC_data %>% filter(serotype %in% vec_serotype) %>% 
    mutate(serotype = fct_relevel(as.factor(serotype), vec_serotype),
           LogGMC = log(GMC),
           LogUCL = log(upper_limit),
           LogLCL = log(lower_limit),
           se_LogGMC = (LogUCL - LogLCL) / (1.96*2))
  return(df)
}
################################################################################
# End of function
################################################################################


# Rearrange the GMC data for vaccine2 & vaccine1 to get GMR
get_GMR <- function(df, vaccine2, vaccine1, vec_st){
  
  # Get ratio for PCV2:PCV1
  df_pcv2 <- subset(df, vaccine == vaccine2) %>% dplyr::select(study_id, serotype, LogGMC, se_LogGMC)
  names(df_pcv2) <- c("study_id","serotype", "LogGMC2", "se_LogGMC2")
  df_pcv1 <- subset(df, vaccine == vaccine1) %>% dplyr::select(study_id, serotype, LogGMC, se_LogGMC)
  names(df_pcv1) <- c("study_id", "serotype", "LogGMC1", "se_LogGMC1")
  
  # Combine into GMR
  df_GMR <- left_join(df_pcv2, df_pcv1, by = c("study_id", "serotype")) %>%
    mutate(serotype = fct_relevel(as.factor(serotype), vec_st),
           LogGMR = LogGMC2 - LogGMC1,
           se_LogGMR = sqrt(se_LogGMC2^2 + se_LogGMC1^2)) %>%
    unique() %>% 
    mutate(uci = LogGMR + 1.96*se_LogGMR,
           lci = LogGMR - 1.96*se_LogGMR)
  return(df_GMR)
}
################################################################################
# End of function
################################################################################


# Create weights
create_wt <- function(df_GMR){
  # Remove -Inf (which will cause NaN), which are from those with mean GMC & lower limit = 0 
  df_GMR <- df_GMR %>% 
    filter(LogGMC1 != -Inf, se_LogGMC1 != -Inf)
  # get inverse-variance weights
  df_wt <- df_GMR %>%
    group_by(study_id, serotype) %>% 
    summarise(variance = se_LogGMR^2,
              inv_var = 1/variance) %>% 
    group_by(serotype) %>%
    mutate(sum_inv_var = sum(inv_var)) %>% 
    ungroup() %>% 
    mutate(wt_inv = inv_var/ sum_inv_var)
  return(df_wt)
}
################################################################################
# End of function
################################################################################


# Get inverse-variance weighted average and its variance
get_wtGMR <- function(df_GMR, df_wt){
  # Remove -Inf (which will cause NaN), which are from those with mean GMC (or lower limit) = 0 
  df_GMR <- df_GMR %>% 
    filter(LogGMC1 != -Inf, se_LogGMC1 != -Inf)
  # Combine weight & GMR data frames and calculate weighted average
  df_wtGMR <- df_GMR %>% left_join(df_wt, by = c("study_id", "serotype")) %>% 
    mutate(wtLogGMR = LogGMR*wt_inv) %>% 
    group_by(serotype) %>%
    summarize(wtAvg = sum(wtLogGMR),
              variance = 1/sum_inv_var,
              se = sqrt(variance)) %>% unique()
  # Rename columns
  names(df_wtGMR) <- c("Serotype", "LogGMR", "var_LogGMR", "se_LogGMR")
  return(df_wtGMR)
}
################################################################################
# End of function
################################################################################


# Patch NA for serotypes with no logGMR info at all once -Inf are removed (can happen when no. of studies included are small)
patch_NA <- function(df_wtGMR, vec_st){
  n_stNA <- sum(vec_st %in% unique(df_wtGMR$Serotype) == FALSE)
  if(n_stNA > 0){
    index_stNA <- which(vec_st %in% unique(df_wtGMR$Serotype) == FALSE)
    vec_stNA <- vec_st[index_stNA]
    df_stNA <- data.frame(Serotype = vec_stNA, LogGMR = NA, var_LogGMR = NA, se_LogGMR = NA)
    df_wtGMR <- rbind(df_wtGMR, df_stNA) %>% 
      mutate(Serotype = as.factor(Serotype),
             Serotype = fct_relevel(Serotype, vec_st)) %>% 
      arrange(Serotype)
    return(df_wtGMR)
  } 
  else{
    print(paste0("n_stNA = ", n_stNA))
    return(df_wtGMR)
  }
}
################################################################################
# End of function
################################################################################


# Combine post-primary & post-boost wtGMR
bind_df <- function(df_postprim, df_postpost){
  df <- rbind(df_postprim %>% mutate(Response = "Post-primary"), 
              df_postpost %>% mutate(Response = "Post-booster")) %>%
    mutate(Response = fct_relevel(as.factor(Response),
                                  c("Post-primary", "Post-booster")),
           uci = LogGMR + 1.96*se_LogGMR,
           lci = LogGMR - 1.96*se_LogGMR)  
}
################################################################################
# End of function
################################################################################
