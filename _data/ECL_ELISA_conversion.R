################################################################################
# To Merck ECL to ELISA, use formula from Ryman et al. 2022 (DOI: 10.1080/14760584.2022.2112179):
# ELISA=(ECL/exp(yintECL))^(1/slopeECL) 
# 
# slopeECL can be found in Table 2 of Nolan et al. 2020 (doi.org/10.4155/bio-2020-0023)
# yintECL can be calculated by:
#
# Substituting 0.35 for ELISA, folddiff*0.35 for ECL, and rearranging terms:
# yintECL=log(folddiff*0.35/(0.35^(slopeECL/1)))
################################################################################
library(pdftools)
library(tidyverse)
library(openxlsx)


# Get slope and intercept from Nolan et al. 2020 Table 2 ------------------
data <- pdf_text("./_data/Nolan_2020_Table2.pdf")

# Extract and tide relevant data section
table <- str_split(data, "\n", simplify = TRUE)
table <- table[c(11:25)]
table <- str_replace_all(table, "\\s{2,}", "|")

# Add col name
col_nm <- "serotype|N|slope|slope_ci|pctdiff|pctdiff_ci|folddiff|agreecoef1|agreecoef2|agreecoef3"
table <- c(col_nm, table[1:15])

# Combine into a table
text_con <- textConnection(table)
data_table <- read.csv(text_con, sep = "|")
data_table <- data_table %>% 
  separate(slope_ci, into = c("slope_lci", "slope_uci"), sep = "\\s*(–|-)\\s*") %>% 
  mutate(intercept = log(folddiff*0.35/(0.35^(slope/1)))) %>% 
  select(serotype, slope, intercept)

# save slope & intercept
saveRDS(data_table, "./_data/ECL_ELISA_conversion.RDS")


# Transform GMC data in PCV15v13 ------------------------------------------

# Set vector of serotypes of interest
vec_7st <- c("4", "6B", "9V", "14", "18C", "19F", "23F")
vec_pcv13.6st <- c("1", "3", "5", "6A", "7F", "19A")
vec_13st <- c(vec_7st, vec_pcv13.6st)

# Get slope and intercept of ECL to ELISA conversion
ECL_ELISA_conversion <- readRDS("./_data/ECL_ELISA_conversion.RDS") %>% 
  filter(serotype %in% vec_13st)

# Get PCV15v13 trial GMC data
# df_15v13 <- read_excel("_data/df_15v13_postprim_n9.xlsx")
df_15v13 <- read_excel("_data/df_15v13_postboost_n11.xlsx")
col_nm <- names(df_15v13)

# Transform PCV15v13 trial GMC data from ECL to ELISA
df_15v13 <- df_15v13 %>% 
  pivot_longer(GMC:lower_limit, names_to = "measurement", values_to = "ECL") %>% 
  left_join(ECL_ELISA_conversion, by = "serotype") %>% 
  mutate(ELISA = round((ECL/exp(intercept))^(1/slope), 2)) %>% 
  select(study_id,participants,vaccine,serotype,measurement, ELISA)

# Rearrange dtat frame and update col names
ls_measurement <- split(df_15v13, df_15v13$measurement)
df_15v13 <- ls_measurement %>% 
  reduce(left_join, by = c("study_id","participants","vaccine","serotype")) %>% 
  select(study_id,participants,vaccine,serotype,ELISA.x, ELISA.y, ELISA)
names(df_15v13) <- col_nm

# Save data
# write.xlsx(df_15v13, "_data/df_15v13_postprim_n9_ELISA.xlsx")
write.xlsx(df_15v13, "_data/df_15v13_postboost_n11_ELISA.xlsx")
