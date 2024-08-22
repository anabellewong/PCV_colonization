################################################################################
# Function for getting colonization observations and censoring subsequent same-
# serotype colonization such that only new acquisitions are included 
################################################################################
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