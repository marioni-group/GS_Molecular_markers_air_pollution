#MAJA covariate preparation

#libraries
library(tidyverse)
library(data.table)
library(VIM)


#1. Load data: methylation IDs, pollution IDs, target file with covariates, simd file
#2. Filter data-set to pollutant/methylation data & check no withdrawals are included
#3. BMI filtering (16 - 50)
#4. Check other covariates 
#5. Log transformation where required
#6. Scale covariates
#7. Impute missing data (KNN)
#8. Save out scaled covariates and save out list of IDs for MAJA. 

#1. Load Data: ids from previous pollutant processing and meth preparation with no covariates
#ids from previous pollutant processing and meth preparation with no covariates for comparison
poll_order <- fread("id_order.csv")
meth_order <- fread("order_resid_scaled_methylation_data.csv")

#load ids for available data
poll_ids <- fread("updated_365_avg_exp_foranalysis_scaleonly.csv")
meth_ids <- readRDS("GS20k_Targets_18869.rds")
covariates <- readRDS("GS_phenos_internal_23Aug2023_REM.rds")
simd <- fread("2024-10-18_simd09_all_categories.csv")

ids_overlap <- intersect(as.character(poll_ids$ID), meth_ids$Sample_Name) 

#check for withdrawals
#Check withdrawals
id_list <- fread("2024-10-31_meth_ids.csv")
ids_overlap2 <- ids_overlap[ids_overlap %in% id_list$V1] 

#2. Filter data
#filter to Ids with methylation and pollution data
cov <- covariates %>%
  dplyr::filter(id %in% ids_overlap2)

#select covariates of interest
cov2 <- cov %>%
  dplyr::select(id, age, sex, bmi, pack_years, units)

#add simd
simd <- simd %>% dplyr::select(id, simd2009v2rank)
cov3 <- left_join(cov2, simd, by = "id")
dim(cov3)

#add batch & SSIDs
batch_target <- meth_ids %>% dplyr::select(Sample_Name, Sample_Sentrix_ID, Batch)
batch_target <- batch_target %>% dplyr::rename(id = Sample_Name)
cov3$id <- as.character(cov3$id)
cov4 <- left_join(cov3, batch_target, by = "id")
dim(cov4) #18556


#3. BMI filtering (16 - 50) 
#filter to bmi range
cov_filt <- cov4 %>%
  dplyr::filter(bmi >=16 & bmi <=50 | is.na(bmi))
dim(cov_filt) #check numbers -> put into analysis document 18512


#4. Check other covariates - any further filtering required?
#check distributions
hist(cov_filt$bmi)
hist(cov_filt$units)
hist(cov_filt$pack_years)
hist(cov_filt$simd2009v2rank)

#look for unrealistic values
cov_filt %>%
  ggplot(aes(y = bmi)) +
  geom_boxplot()

cov_filt %>%
  ggplot(aes(y = units)) +
  geom_boxplot()

#write out with no imputation
fwrite(cov_filt, file = "cov18512_noimpt.csv")

#5. Log transformation: log transform bmi, pack_years, units

covf <- cov_filt %>%
  mutate(logbmi = log(bmi),
         logunits = log(units + 1),
         logsmok = log(pack_years + 1))
sapply(covf, function(y) sum(length(which(is.na(y)))))

#6. Impute missing data (KNN)
#Prep for scaling and knn imputation
 #Female is 1
covf <- covf %>%
  mutate(sex = ifelse(sex == "F", 1, 0))

#select transformed variables only prior to imputation
covf <- covf %>% dplyr::select(id, Sample_Sentrix_ID, age, sex, Batch, logsmok, logbmi, logunits, simd2009v2rank)
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
fwrite(covf4, file = "cov18512_sc_impt.csv")

#save out covariates file minus smoking for use in methylation data pre-prep
covf5 <- covf4 %>% dplyr::select(-logsmok)
fwrite(covf5, file ="cov18512_nosmoksc_impt.csv")

