#Multivariate analysis of proteins and pollutants

library(tidyverse)
library(data.table)
library(here)
library(VIM)

#1. Check overlapping Ids of protein & pollution data
#2. Filter to covariates to these IDs
#3. Filter covariates by bmi
#4. Impute missing covariates
#5. Regress protein data
#6. Run Maja

#1. Find IDs with protein data and pollution data
proteins <- readRDS(here("GS_ProteinGroups_RankTransformed_23Aug2023.rds"))
pollution <- fread("updated_365_avg_exp_foranalysis_scaleonly.csv")
dim(proteins) 
dim(pollution) 

length(proteins$id %in% pollution$ID) 
proteins2 <- proteins %>%
 dplyr::filter(id %in% pollution$ID) 

 #check withdrawals
non_drop_outs <- read.table("GS_ids_20241007.txt")
proteins3 <- proteins2 %>%
  dplyr::filter(id %in% non_drop_outs$V1) 

#Filter covariates to these Ids
covariates <- readRDS("GS_phenos_internal_23Aug2023_REM.rds")  
cov <- covariates %>%
  dplyr::filter(id %in% proteins3$id) 

#select covariates of interest
cov2 <- cov %>%
  dplyr::select(id, age, sex, bmi, pack_years, units)

#add simd
simd <- fread("2024-10-18_simd09_all_categories.csv")
simd <- simd %>% dplyr::select(id, simd2009v2rank)
cov3 <- left_join(cov2, simd, by = "id")
dim(cov3)

#add batch & SSIDs
meth_ids <- readRDS("GS20k_Targets_18869.rds")
batch_target <- meth_ids %>% dplyr::select(Sample_Name, Sample_Sentrix_ID)
batch_target <- batch_target %>% dplyr::rename(id = Sample_Name)
cov3$id <- as.character(cov3$id)
cov4 <- left_join(cov3, batch_target, by = "id")
dim(cov4) 


#3. BMI filtering (16 - 50)
#filter to bmi range
cov_filt <- cov4 %>%
  dplyr::filter(bmi >=16 & bmi <=50 | is.na(bmi))

sapply(cov_filt, function(y) sum(length(which(is.na(y)))))

#write out with no imputation
fwrite(cov_filt, file = "proteins_pollution/cov_noimp") 

#5. Log transformation: log transform bmi, pack_years, units

covf <- cov_filt %>%
  mutate(logbmi = log(bmi),
         logunits = log(units + 1),
         logsmok = log(pack_years + 1))
sapply(covf, function(y) sum(length(which(is.na(y)))))

covf <- covf %>% dplyr::select(-Sample_Sentrix_ID)

#6. Impute missing data (KNN)
#Prep for scaling and knn imputation
 #Female is 1
covf <- covf %>%
  mutate(sex = ifelse(sex == "F", 1, 0))

#select transformed variables only prior to imputation
covf <- covf %>% dplyr::select(id, age, sex, logsmok, logbmi, logunits, simd2009v2rank)
str(covf)
covf$sex <- as.factor(covf$sex)

##6. Scale covariates: Scale numeric covariates
covf2 <- covf %>% modify_if(., is.numeric, scale)
sapply(covf2, function(y) sum(length(which(is.na(y)))))

#7. Impute missing data (KNN)
set.seed(1.234)
covf3 <- kNN(covf2, k = 10)
sapply(covf3, function(y) sum(length(which(is.na(y))))) #0

#remove imputation variables created
covf4 <- subset(covf3, select = id:simd2009v2rank)

#Save out covariates file for reference.
fwrite(covf4, file = "cov15314_sc_impt.csv")
#save out covariates file minus smoking for use in methylation data pre-prep
covf5 <- covf4 %>% dplyr::select(-logsmok)
fwrite(covf5, file ="cov15314_nosmoksc_impt.csv")
