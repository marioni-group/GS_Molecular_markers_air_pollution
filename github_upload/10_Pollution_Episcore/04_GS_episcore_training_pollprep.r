#Pollution data-prep

library(tidyverse)
library(data.table)

#Load un-transformed data
pollution <- fread("annual_avg_exp_foranalysis_ntr.csv")

#Transpose
pollution <- pollution %>%
  dplyr::select(ID, Pollutant, days_365_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )
dim(pollution)
# [1] 22071     9

#Subset to IDs in prepped target_filt file (EPIC_poll_episcore_meth_dataprep)
tfilt <- fread("/target_15516_pollution_episcore_input.csv")
# > dim(tfilt)
# [1] 15516    27

poll <- pollution %>% filter(ID %in% tfilt$id)
dim(poll)
# [1] 15516     9

#Check orders match tfilt
poll <- poll[match(tfilt$id, poll$ID),]
identical(poll$id, tfilt$id) #FALSE
poll <- poll %>% rename(id = ID)
poll$id <- as.character(poll$id)
tfilt$id <- as.character(tfilt$id)
identical(poll$id, tfilt$id) #TRUE

#add smoking
covariates <- fread("cov18512_noimpt.csv")
cov <- covariates %>% filter(id %in% tfilt$id)
tfilt$id <- as.character(tfilt$id)
cov$id <- as.character(cov$id)
cov <- cov[match(tfilt$id, cov$id),]
identical(cov$id, tfilt$id) #TRUE
sum(is.na(cov$pack_years)) 
cov <- cov %>% mutate(logsmok = log(pack_years + 1)) %>% select(-pack_years)

#scale before imputation
str(cov)
cov <- cov %>% modify_if(is.numeric, scale)
cov2 <- cov

#remove ID column
cov2 <- cov2 %>% select(-Sample_Sentrix_ID)
cov2 <- as.data.frame(cov2)

#Impute smoking based on other covariates
cov2 <- kNN(cov2, variable = "logsmok", k = 10)
smok <- cov2 %>% select(id, logsmok)
sum(is.na(smok))

#join 
poll <- poll %>% modify_if(is.numeric, scale)
poll <- left_join(poll, smok, by = "id")
dim(poll)
head(poll)

sum(is.na(poll))

#for comparison
smok <- fread("cov18512_sc_impt.csv")
smok <- smok %>% filter(id %in% poll$id)
identical(smok$id, poll$id)
smok <- smok[match(poll$id, smok$id),]
smok$id <- as.character(smok$id)
identical(smok$id, poll$id) #TRUE
cor.test(smok$logsmok, poll$logsmok)
# sample estimates:
#       cor 
# 0.9967616 

#check structure
str(poll)

#Save
fwrite(poll, file = "poll_scaled_15516.csv")

#save for maja
poll2 <- poll %>% select(-id)
fwrite(poll2, file = "poll_scaled_15516.txt", sep = " ", col.names = FALSE)

#Double check orders against methylation data prep
cov <- fread("target_15516_pollution_episcore_input.csv")
cov$id <- as.character(cov$id)
identical(cov$id, poll$id) #TRUE
