#3 pollution data prep for protein analysis
library(data.table)
library(tidyverse)

#1. Filter to correct IDs
pollution <- fread("annual_avg_exp_foranalysis_scaled.csv")

cov <- fread("cov15314_sc_impt.csv")
pollution <- pollution %>% 
  dplyr::filter(ID %in% cov$id)

#check order
identical(pollution$ID, cov$id) #FALSE
pollution <- pollution %>%
  dplyr::rename(id = ID)

#Match orders
pollution <- pollution[match(cov$id, pollution$id),]  
identical(pollution$id, cov$id) #TRUE

#add smoking
smok <- cov %>% dplyr::select(id, logsmok)
poll <- left_join(pollution, smok, by = "id")

#save out with ids and colnames 
fwrite(poll, file = "poll15314_sc.csv" )


#save out for use with maja
poll <- poll %>%
  dplyr::select(-c(id))

fwrite(poll, file = "poll15314_sc.txt",  sep = " ", col.names = FALSE)
