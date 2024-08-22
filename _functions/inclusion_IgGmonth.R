################################################################################
# Function for including subjects based on exposure measurement:
# whether IgG measurement at the selected month (7m or 13m) is available
################################################################################
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