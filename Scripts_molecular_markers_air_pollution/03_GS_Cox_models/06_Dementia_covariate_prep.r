#Covariate prep: Dementia


#For all models
#1. Age + sex + kinship as covariates
#2. Age + sex + smoking + bmi + alcohol units + SIMD rank + kinship
#3. Age + sex + smoking + bmi + alcohol units + SIMD rank + APOE status + kinship (Coxme model)


#_____________________________________________________________________________
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
pheno <- readRDS("GS_phenos_internal_23Aug2023_REM.rds") 
cox <- fread("dementia_coxdf_10761.csv") 
simd <- fread("2024-10-18_simd09_all_categories.csv") 
target <- readRDS("GS20k_Targets_18869.rds") 
target_ids <- target %>% dplyr::select(Sample_Name, Sample_Sentrix_ID)
include_ids <- cox %>% dplyr::select(id)
include_ids <- include_ids %>% dplyr::rename(Sample_Name = id)
include_ids$Sample_Name <- as.character(include_ids$Sample_Name)
pollution <- fread("annual_avg_exp_foranalysis_ntr.csv")
diseases <- fread("2025-01-09_diseases.csv")

#add SSIDs to IDs for filtering/order matching
both_ids <- left_join(include_ids, target_ids, by = "Sample_Name")
dim(both_ids)
identical(both_ids$Sample_Name, include_ids$Sample_Name) #TRUE: order remains the same

#Filter to Ids to include
pheno_filtered <- pheno %>%
  dplyr::filter(id %in% include_ids$Sample_Name) 

#now filter by bmi
pheno_bmifilt <- pheno_filtered %>% dplyr::filter(bmi >= 16 & bmi <=50 | is.na(bmi)) 
summary(pheno_bmifilt$bmi)

#filter coxdf by IDs - get DF with final numbers for analysis
cox2 <- cox %>% dplyr::filter(id %in% pheno_bmifilt$id) 

fwrite(cox2, file = "dementia_coxdf_140525_10739.csv")

#___________________________________________________________
#create df of phenotypes
#3. Age + sex + smoking + bmi + alcohol units + SIMD rank + APOE status + kinship (Coxme model) + diabetes + hypertension + 

cox_phenos <- left_join(pheno_bmifilt, simd, by = "id")

#Create apoe variable 
cox_phenos <- cox_phenos %>%
  mutate(e4_count = case_when(
    apoe == "e4e4" ~ 2,
    str_detect(apoe, "e4") ~ 1,
    TRUE ~ 0))

table(cox_phenos$e4_count)

cox_phenos <- cox_phenos %>% dplyr::select(id, age, sex, y_appt, pack_years, 
                                            bmi, units, e4_count, simd2009v2rank, 
                                            HDL_cholesterol, Total_cholesterol)
str(cox_phenos)

dis <- diseases %>% 
   filter(disease == "diabetes" | disease == "hypertension") %>%
   filter(incident == "0" & source == "Secondary_Care") %>%
   rename(dt_diag = dt1_ym)

dis <- dis %>%
  mutate(diabetes = ifelse(disease == "diabetes" & incident == 0, 1, 0),
         hypertension = ifelse(disease == "hypertension" & incident == 0, 1, 0))
  
  
dis2 <- dis %>%
 group_by(id) %>%
 summarise(
  diabetes = max(as.numeric(diabetes)),
  hypertension = max(as.numeric(hypertension))
 ) %>% ungroup()

#add to cox_phenos
cox_phenos2 <- left_join(cox_phenos, dis2, by = "id")
cox_phenos2 <- cox_phenos2 %>%
  mutate(diabetes = coalesce(diabetes, 0),
        hypertension = coalesce(hypertension, 0))

#save out
fwrite(cox_phenos2, file = "Dementia/pheno_bmifilt_080725.csv")

#Create summary table
#Extract demographic info from non-imputed covariates
#Add cox data

cox_phenos2 <- cox_phenos2 %>% dplyr::select(-c(age, sex, y_appt))
cox_df <- left_join(cox2, cox_phenos2, by = "id")


#add pollutants too
pollution2 <- pollution %>% 
  dplyr::rename(id = ID) %>%
  dplyr::filter(id %in% cox_df$id) %>%
  dplyr::select(id, Pollutant, days_365_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )

cox_df2 <- left_join(cox_df, pollution2, by = "id") 

#Create summary table
library(arsenal)
table <- tableby(dementia ~ sex + tte_years + is_death + age + bmi + pack_years + units + simd2009v2rank + 
                as.factor(e4_count) + as.factor(y_appt) + HDL_cholesterol + Total_cholesterol + 
                as.factor(diabetes) + as.factor(hypertension) +
                PM25 + NO2 + NO + O3 + PM10 + SO2 + NO3_C + NO3_F, 
                  data = cox_df2, 
                  test = FALSE, 
                  control = tableby.control(
                  numeric.stats = c("mean", "sd", "median", "q1q3", "Nmiss2")))
tab1 <- as.data.frame(summary(table))

fwrite(tab1, file = "Dementia/summary_demo_andpoll_10739_080725.csv") #bmi filtering, no imputation


#Impute
#Prep for scaling and knn imputation
str(cox_df2)
sapply(cox_df2, function(y) sum(length(which(is.na(y))))) #checking n missing

#drop cox outcome variables for imputation
cox_outcomes <- cox_df2 %>% dplyr::select(id, dementia, is_death, tte_years, status)
covariates <- cox_df2 %>% dplyr::select(id, age, sex, y_appt, pack_years, bmi, units, e4_count, simd2009v2rank, HDL_cholesterol, Total_cholesterol, diabetes, hypertension)

#Set up data
covariates$id <- as.character(covariates$id)
covariates$sex <- as.factor(covariates$sex)
covariates$e4_count <- as.factor(covariates$e4_count)
covariates$y_appt <- as.factor(covariates$y_appt)
covariates$logsmok <- log10(covariates$pack_years + 1)
covariates$logunits <- log10(covariates$units + 1)
covariates$logbmi <- log10(covariates$bmi)
covariates$age <- as.character(covariates$age) #set as character so can centre without scaling
covariates$diabetes <- as.factor(covariates$diabetes)
covariates$hypertension <- as.factor(covariates$hypertension)
str(covariates)

#Scale - required prior to knn imputation
cov2 <- covariates %>% modify_if(., is.numeric, scale)
cov2$age <- scale(as.numeric(cov2$age), scale = FALSE)
str(cov2)

#Check if any rows/columns have missing data > 0.4
cov3 <- cov2 %>% dplyr::select(-c(pack_years, units, bmi))
max_missing <- 0.6
rows_with_high_missing <- which(rowMeans(is.na(cov3)) > max_missing) 
rows_with_high_missing 
cols_with_high_missing <- which(colMeans(is.na(cov3)) > max_missing) #0
cols_with_high_missing

#Impute missing data
library(VIM)
set.seed(12)
#drop non-logged data

cov_imp <- VIM::kNN(as.data.frame(cov3), k = 5)
col_tokeep <- colnames(cov_imp)[-1]

#Rejoin with prepared cox_df 
cox_df3 <- cox_df2 %>% dplyr::select(-any_of(c("pack_years", "age", "sex", "e4_count", 
                                                "units", "bmi", "simd2009v2rank", "HDL_cholesterol", 
                                                "Total_cholesterol", "diabetes", "hypertension"))) #remove columns to be replaced with imputed
cox_df3$id <- as.character(cox_df3$id)
cox_df3 <- left_join(cox_df3, cov_imp, by = "id") #add in imputated data
str(cox_df3)
cox_df3 <- cox_df3 %>% select(-matches("_imp"))
cox_df3 <- cox_df3 %>% select(-y_appt.y) %>% rename(y_appt = y_appt.x)
str(cox_df3)


pollutants <- unique(pollution$Pollutant)#list of pollutants
cox_df4 <- cox_df3 %>% 
  mutate(across(all_of(pollutants), ~scale(.)))
str(cox_df4)

#save out for use
fwrite(cox_df4, file = "Dementia/cox_df_sc_imp_08072025.csv")

