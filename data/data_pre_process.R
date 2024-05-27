rm(list = ls())
##------------------------------------------------------------------------------
## Load libraries 
##------------------------------------------------------------------------------
library(haven)

##------------------------------------------------------------------------------
## Data pre-processing
##------------------------------------------------------------------------------
## Download original data from 
## https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN1/22633&version=1.1 
## Save 090113_TotWatDat_cor_merge_Price.dta file in data folder

original_data = read_dta("090113_TotWatDat_cor_merge_Price.dta")

original_data %>% mutate(Y = jun07_3x+jul07_3x+aug07_3x+sep07_3x, D = treatment) %>%
  select(Y,D,jun06,jul06,aug06,sep06,oct06,nov06,dec06,jan07,feb07,mar07,apr07_3x,may07_3x) %>%
  write_csv("data_ferraroprice.csv")

