#CVD cox_models covariate prep

#For all models - run as coxme
      # paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      # paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + logunits + logbmi + logsmok + simd2009v2rank"), 
      # paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id) + logunits + logbmi + logsmok + simd2009v2rank + as.factor(diabetes_Y) + Total_cholesterol + as.factor(high_BP_Y)"")
 

#Following exploratory analysis, will knn impute missing covariates to avoid data loss and retain power
#Filter by age 40-69
#Filter out BMI extremes


#Load Packages
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("impute")
library(impute)
library(here)
library(tidyverse)
library(data.table)
library(VIM)

#________________________________________________________________________________
#Load files
pheno <- readRDS("/GS_phenos_internal_23Aug2023_REM.rds") 
simd <- fread("2024-10-18_simd09_all_categories.csv") 
target <- readRDS("GS20k_Targets_18869.rds") #SSIDs and IDs
target_ids <- target %>% dplyr::select(Sample_Name, Sample_Sentrix_ID)
diseases <- fread("2025-01-09_diseases.csv")
pollution <- fread("annual_avg_exp_foranalysis_ntr.csv")

#Load cox dataframe of interest 

filename <- "CVD_htn_basic_coxdf_agebmifilt.csv" #Generated in 01_CVD_coxmodeldataprep.r
path <-  ""
cox <- fread(paste0(path,filename))

#Generate list of IDs to include for phenotypes
include_ids <- cox %>% dplyr::select(id)
include_ids <- include_ids %>% dplyr::rename(Sample_Name = id)
include_ids$Sample_Name <- as.character(include_ids$Sample_Name)

#add SSIDs to IDs for filtering/order matching
both_ids <- left_join(include_ids, target_ids, by = "Sample_Name")
dim(both_ids)
identical(both_ids$Sample_Name, include_ids$Sample_Name) #TRUE: order remains the same

#Filter to Ids to include
pheno_filtered <- pheno %>%
  dplyr::filter(id %in% include_ids$Sample_Name) 
summary(pheno_filtered$bmi)

#___________________________________________________________
#create df of phenotypes

cox_phenos <- left_join(pheno_filtered, simd, by = "id")
cox_phenos <- cox_phenos %>% dplyr::select(id, age, sex, y_appt, pack_years, bmi, units, simd2009v2rank, HDL_cholesterol, Total_cholesterol)
str(cox_phenos)

#add in htn and diabetes flags
coxdf <- fread("CVD_htn_basic_coxdf_agebmifilt.csv")
add <- coxdf %>% select(id, diabetes, hypertension)

cox_phenos <- left_join(cox_phenos, add, by = "id")

#Create summary table
#Extract demographic info from non-imputed covariates
#add cox data

cox_cvd <- coxdf %>% dplyr::select(id, cvd, status_IHD, status_MI,
                                    status_stroke, cardiac_death, dead, 
                                    composite_cardiac_outcome, tte_years, 
                                    status) %>%
  dplyr::filter(id %in% cox_phenos$id)
cox_cvd2 <- left_join(cox_cvd, cox_phenos, by = "id")


#add pollutants too
pollution2 <- pollution %>% 
  dplyr::rename(id = ID) %>%
  dplyr::filter(id %in% cox_cvd2$id) %>%
  dplyr::select(id, Pollutant, days_365_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )


cox_cvd2 <- left_join(cox_cvd2, pollution2, by = "id") 

#Create summary table
library(arsenal)
table <- tableby(composite_cardiac_outcome ~ sex + tte_years + as.factor(cvd) + as.factor(cardiac_death) + as.factor(dead) + as.factor(status_IHD) + as.factor(status_MI) + as.factor(status_stroke) + age + bmi + pack_years + units + simd2009v2rank + as.factor(y_appt) + PM25 + NO2 + NO + O3 + PM10 + SO2 + NO3_C + NO3_F, 
                  data = cox_cvd2, 
                  test = FALSE, 
                  control = tableby.control(
                  numeric.stats = c("mean", "sd", "median", "q1q3", "Nmiss2")))
tab1 <- as.data.frame(summary(table))

fwrite(tab1, file = "summary_cvd_andpoll_13023_070725.csv") # age & bmi filtering, no imputation



#__________________________________________________________________________________________
#Impute

#Prep for scaling and knn imputation
str(cox_df2)
sapply(cox_df2, function(y) sum(length(which(is.na(y))))) #checking n missing - below is for full dataset

#set up for imputation
covariates <- cox_phenos
covariates$id <- as.character(covariates$id)
covariates$sex <- as.factor(covariates$sex)
covariates$y_appt <- as.factor(covariates$y_appt)
covariates$logsmok <- log10(covariates$pack_years + 1)
covariates$logunits <- log10(covariates$units + 1)
covariates$logbmi <- log10(covariates$bmi)
covariates$age <- as.character(covariates$age) #so can centre without scaling
covariates$diabetes <- as.character(covariates$diabetes)
covariates$hypertension <- as.character(covariates$hypertension)
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
cols_with_high_missing <- which(colMeans(is.na(cov3)) > max_missing)
cols_with_high_missing

#Impute missing data
library(VIM)
set.seed(12)

variables <- c("simd2009v2rank", "logsmok", "logunits", "logbmi", "HDL_cholesterol", "Total_cholesterol")
cov_imp <- VIM::kNN(as.data.frame(cov3), variable = variables,  k = 5)
cox_cvd$id <- as.character(cox_cvd$id)
cov_imp <- cov_imp %>% select(-matches("_imp"))
cox_cvd3 <- left_join(cox_cvd, cov_imp, by = "id")
str(cox_cvd3)

#add in pollutants
pollutants <- colnames(pollution2)[-1]#list of pollutants
pollution2$id <- as.character(pollution2$id)
cox_cvd4 <- left_join(cox_cvd3, pollution2, by = "id")
cox_cvd4 <- cox_cvd4 %>% mutate(across(all_of(pollutants), ~scale(.)))
str(cox_cvd4)


#add in extra outcomes prepped in 01_coxph_dataprep_CVD

agebmifilt <- fread("CVVD_basic_coxdf_agebmifilt.csv")
outcomes <- agebmifilt %>% select(id, status_IHD, status_MI, status_stroke)

coxdf <- fread("cvd_cox_df_agebmifilt_sc_imp.csv")
coxdf <- left_join(coxdf, outcomes, by = "id")

#set NAs to 0 
coxdf2 <- coxdf %>%
  mutate(across(c(status_IHD, status_MI, status_stroke), ~ replace_na(., 0)))

#save out
fwrite(coxdf2, "cvd_cox_df_agebmifilt_sc_imp_extraoutcomes.csv")


