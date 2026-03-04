library(impute)
library(here)
library(tidyverse)
library(data.table)
library(VIM)

#_______________________________________________________________________________________________________________________________________________
#Covariates for HTN coxdf 
#_______________________________________________________________________________________________________________________________________________
pheno <- readRDS("GS_phenos_internal_23Aug2023_REM.rds") 
target <- readRDS("GS20k_Targets_18869.rds") #SSIDs and IDs
simd <- fread("2024-10-18_simd09_all_categories.csv") 
target_ids <- target %>% dplyr::select(Sample_Name, Sample_Sentrix_ID)

filename <- "CVD_htn_basic_coxdf_bmifilt.csv"
path <-  ""

cox <- fread(paste0(path,filename))
dim(cox)
cox_htn <- cox %>% filter(!hypertension == 1)
cox %>% filter(hypertension == 1) %>% dim() 
dim(cox_htn) 


include_ids <- cox_htn %>% dplyr::select(id)
include_ids <- include_ids %>% dplyr::rename(Sample_Name = id)
include_ids$Sample_Name <- as.character(include_ids$Sample_Name)
pollution <- fread("annual_avg_exp_foranalysis_ntr.csv")

#add SSIDs to IDs for filtering/order matching
both_ids <- left_join(include_ids, target_ids, by = "Sample_Name")
dim(both_ids) 
identical(both_ids$Sample_Name, include_ids$Sample_Name) #TRUE: order remains the same

#Filter to Ids to include
pheno_filtered <- pheno %>%
  dplyr::filter(id %in% include_ids$Sample_Name) 
summary(pheno_filtered$bmi) #already filtered by BMI


#___________________________________________________________
#create df of phenotypes
#3. Age + sex + smoking + bmi + alcohol units + SIMD rank + y_appt + kinship (Coxme model)

cox_phenos <- left_join(pheno_filtered, simd, by = "id")

cox_phenos <- cox_phenos %>% dplyr::select(id, age, sex, y_appt, pack_years, bmi, units, simd2009v2rank, HDL_cholesterol, Total_cholesterol)
str(cox_phenos)

#save out
fwrite(cox_phenos, file = "pheno_bmifilt_forhtn_070725.csv") 

#add in htn and diabetes flags
add <- cox_htn %>% select(id, diabetes, hypertension)
cox_phenos <- left_join(cox_phenos, add, by = "id")

#Create summary table
#Extract demographic info from non-imputed covariates
#Add cox_data

cox_htn2 <- cox_htn %>% dplyr::select(id, htn, status_htn, dead, composite_cardiac_outcome, tte_htn_years) %>%
  dplyr::filter(id %in% cox_phenos$id)
cox_htn2 <- left_join(cox_htn2, cox_phenos, by = "id")
cox_htn2 <- cox_htn2 %>% select(-hypertension)

#add pollutants too
pollution2 <- pollution %>% 
  dplyr::rename(id = ID) %>%
  dplyr::filter(id %in% cox_htn2$id) %>%
  dplyr::select(id, Pollutant, days_365_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )

cox_htn2$id <- as.character(cox_htn2$id)
pollution2$id <- as.character(pollution2$id)
cox_htn2 <- left_join(cox_htn2, pollution2, by = "id")

#Create summary table
library(arsenal)
table <- tableby(status_htn ~ sex + tte_htn_years + as.factor(dead) + age + bmi + pack_years + units + simd2009v2rank + as.factor(y_appt) + HDL_cholesterol + Total_cholesterol + as.factor(diabetes) + PM25 + NO2 + NO + O3 + PM10 + SO2 + NO3_C + NO3_F, 
                  data = cox_cvd2, 
                  test = FALSE, 
                  control = tableby.control(
                  numeric.stats = c("mean", "sd", "median", "q1q3", "Nmiss2")))
tab1 <- as.data.frame(summary(table))

fwrite(tab1, file = "summary_htn_andpoll_20753_070725.csv") #bmi filtering, no imputation


#Prep for scaling and imputation
covariates <- cox_phenos %>% select(-hypertension)
covariates$id <- as.character(covariates$id)
covariates$sex <- as.factor(covariates$sex)
covariates$y_appt <- as.factor(covariates$y_appt)
covariates$logsmok <- log10(covariates$pack_years + 1)
covariates$logunits <- log10(covariates$units + 1)
covariates$logbmi <- log10(covariates$bmi)
covariates$age <- as.character(covariates$age) #so can centre without scaling
covariates$diabetes <- as.character(covariates$diabetes)
str(covariates)

#Scale - required prior to knn imputation
cov2 <- covariates %>% modify_if(., is.numeric, scale)
cov2$age <- scale(as.numeric(cov2$age), scale = FALSE) #centre but don't scale age
str(cov2)

#Check if any rows/columns have missing data > 0.4
cov3 <- cov2 %>% dplyr::select(-c(pack_years, units, bmi)) #remove phenotypes which have been log transformed
max_missing <- 0.6
rows_with_high_missing <- which(rowMeans(is.na(cov3)) > max_missing) 
rows_with_high_missing 
cols_with_high_missing <- which(colMeans(is.na(cov3)) > max_missing) #0
cols_with_high_missing

#Impute missing data
library(VIM)
set.seed(12)

#impute
variables <- c("simd2009v2rank", "logsmok", "logunits", "logbmi", "HDL_cholesterol", "Total_cholesterol")
cov_imp <- VIM::kNN(as.data.frame(cov3), variable = variables,  k = 5)

#Tidy up
cox_htn2$id <- as.character(cox_htn2$id)
cov_imp <- cov_imp %>% select(-matches("_imp"))
cox_htn2 <- cox_htn %>% select(-c(hypertension, age, sex, diabetes))
cox_htn2$id <- as.character(cox_htn2$id)
cox_htn3 <- left_join(cox_htn2, cov_imp, by = "id") #add in imputated data
str(cox_htn3)

#add in pollutants
pollutants <- colnames(pollution2)[-1]#list of pollutants
pollution2$id <- as.character(pollution2$id)
cox_htn4 <- left_join(cox_htn3, pollution2, by = "id")
cox_htn4 <- cox_htn4 %>% mutate(across(all_of(pollutants), ~scale(.)))
str(cox_htn4)

#save out
fwrite(cox_htn4, file = "htn_cox_df_bmifilt_sc_imp_070725.csv")