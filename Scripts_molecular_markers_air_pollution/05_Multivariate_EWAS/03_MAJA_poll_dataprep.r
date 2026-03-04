#MAJA pollution/outcome data_prep
#For troubleshooting steps/alternative phenotype formats see scripts/pollution_analysis/maja/EWAS/full_run_2/03_MAJA_poll_dataprep_full2812_majatests.r


library(tidyverse)
library(data.table)

#phenotypes = 8 pollutants + smoking
#covariates = age, sex, Batch, logbmi, logunits, simdrank - regressed out of methylation data
#filtered to bmi

#The phenotypic data is required to be stored as txt format with one trait per column, standardized per column. Not available measures need to be imputed.

#Load covariates (pollution)
annual_dft <- fread("annual_avg_exp_foranalysis_scaled.csv")

#add scaled smoking data
smok <- fread("cov18512_sc_impt.csv")
cov <- fread("cov18512_sc_impt.csv")
dim(smok)
smok <- smok %>% dplyr::select(id, Sample_Sentrix_ID, logsmok)

#filter to correct ids
annual_dft <- annual_dft %>% dplyr::filter(ID %in% smok$id)
annual_dft <- annual_dft %>%
  dplyr::rename(id = ID)


#Check withdrawals
id_list <- fread("2024-10-31_meth_ids.csv")
non_drop_outs <- read.table("GS_ids_20241007.txt")
check <- annual_dft %>% dplyr::filter(id %in% non_drop_outs$V1)  

#Join smoking phenotype
annual_dft2 <- left_join(annual_dft, smok, by = "id")


#check orders
order <- fread("SSID_order_mvalinput.csv")
annual_dft2$id <- as.character(annual_dft2$id)
identical(annual_dft2$Sample_Sentrix_ID, order$"result$colnames") #

order_ids <- order %>%
  rename(Sample_Sentrix_ID = "result$colnames")
annual_dft2 <- annual_dft2[match(order_ids$Sample_Sentrix_ID, annual_dft2$Sample_Sentrix_ID), ]
identical(annual_dft2$Sample_Sentrix_ID, order_ids$Sample_Sentrix_ID) 



#Save out with order of pollutants as previous
new_ord2 <- annual_dft2 %>% dplyr::select(id, Sample_Sentrix_ID,PM25,O3,NO,NO2,NO3_C,NO3_F,PM10,SO2,logsmok)
fwrite(new_ord2, file = "full_poll_ids_newphenorder_2905.csv")

#and save out as MAJA input
new_ord3 <- new_ord2 %>% select(-c(id,Sample_Sentrix_ID))
fwrite(new_ord3, file = "full_poll_ids_newphenorder_2905.txt", sep = " ", col.names = FALSE)


#check trait orders
pollids <- fread("full_poll_ids.csv")
input_file <- fread("update_0425_fullpoll.txt")
trait_list_checked <- colnames(pollids)
trait_list_checked <- trait_list_checked[-1]
trait_list_checked <- trait_list_checked[-9]

#check Ids
ssid_order <- fread("SSID_order_mvalinput.csv")
ssid_order <- ssid_order %>% dplyr::rename(SSID = "result$colnames")
ids <- cov %>% dplyr::select(id, Sample_Sentrix_ID)
identical(ssid_order$SSID, ids$Sample_Sentrix_ID) 
identical(pollids$id, ids$id) 
