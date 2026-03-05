#Pollution - CVD
#Code from: /Cluster_Filespace/Marioni_Group/Josie/Proteins/scripts/pollution_analysis/cox_models/CVD/01_coxph_dataprep_CVD_basicdf.r
#Simplified to contain only the data-prep used for the included analysis

#. 1. Set-up incident diseases of interest - generate yes(1) or no(0) column (CVD: diagnosis or diagnosis at death)
#. 2. Remove prevalent disease (individuals diagnosed before onset of study)
#. 3. Filter to non-drop-out IDs
#. 4. Set-up CVD outcomes
#. 5. Set-up cox model variables: event date, TTE, Status
#. 6. Filter to 40 - 69 (SCORE2 age range) - compare numbers pre and post filtering


library(data.table)
library(tidyverse)
library(stringr)
library(here)
library(magrittr)
library(lubridate)
library(here)

#Cox model libraries
library(survival)
library(survminer)

#for plot ordering
library(forcats)

#Load Data:
#Required:
  # Pollution data for IDs
  # Non-drop-outs for IDs
  # Health outcome data - including ICD codes
  # Phenotype data: full covariates to include - 
  # SIMD
  # Deaths data (including cause)
  # GS appointment dates

#Load Data

#Pollution data: 365days for now
annual_df <- fread("/annual_avg_exp_foranalysis_ntr.csv") 
non_drop_outs <- read.table("GS_ids_20241007.txt")
deaths1 <- fread("2025-01-09_deaths_pollids.csv")
pheno <- readRDS("GS_phenos_internal_23Aug2023_REM.rds")
diseases <- fread("2025-01-09_diseases.csv")
gsappt <- fread("GS_appt.txt") 




#________________________________________________________________________________________
#1. Generate outcome df for individuals in test set: cvd cases or not

#Filter pheno to ids with pollution data
pheno$id <- as.character(pheno$id)
pheno <- pheno %>% dplyr::filter(id %in% annual_df$ID) 
check <- annual_df %>% dplyr::filter(ID %in% pheno$id) 

#subset only to columns of interest for creating the cox df
pheno_cox <- pheno %>% dplyr::select(id, age, sex, dob_ym, y_appt)

#______________________________________________
#Set-up CVD outcomes
#subset disease phenotypes to ids in pheno
diseases <- diseases %>% dplyr::filter(.,id %in% pheno$id)
dim(diseases)

#Set-up CVD outcomes
cvd <- c("CHD_NOS", "Isch_stroke", "myocardial_infarction")
cvd_cases <- diseases %>% 
  dplyr::filter(disease %in% cvd & source == "Secondary_Care" & incident == "1")
length(unique(cvd_cases$id)) 

#Look at prevalent disease
prevalent_cvd <- diseases %>%
  dplyr::filter(disease %in% cvd & source == "Secondary_Care" & incident == "0" ) 
length(unique(prevalent_cvd$id)) 

#Exclude prevalent cases and create marker for cvd_diagnosis
cvd_cases <- cvd_cases %>%
  dplyr::filter(!id %in% prevalent_cvd$id) %>%
  mutate(cvd = c("1")) 
length(unique(cvd_cases$id)) 

#also make flags for individual conditions: 
unique(cvd_cases$disease)

cvd_cases <- cvd_cases %>%
  mutate(MI = ifelse(disease == "myocardial_infarction", "1", "0"),
         Isch_stroke = ifelse(disease == "Isch_stroke", "1", "0"),
         CHD_NOS = ifelse(disease == "CHD_NOS", "1", "0"))
table(cvd_cases$CHD_NOS) 
table(cvd_cases$MI) 
table(cvd_cases$Isch_stroke) 
table(cvd_cases$cvd) 


#generate dates and select earliest diagnosis date where multiple diagnoses are present
cvd_cases2 <- cvd_cases %>%
  dplyr::mutate(dt1_ym = ym(dt1_ym),
                gs_appt = ym(gs_appt)) 

cvd_cases3 <- cvd_cases2 %>%
  dplyr::group_by(id) %>%
  dplyr::filter(dt1_ym == min(dt1_ym)) %>% #keep first entry where diagnoses on different dates
  slice_min(order_by = dt1_ym, with_ties = FALSE) %>% #keep first entry where diagnoses on same date
  dplyr::ungroup() 

#generate cvd deaths info & check if any overlap with cvd diagnoses
heart_deaths <- c("I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I11",
         "I13", "I20", "I21", "I22", "I23", "I24", "I25", "I26", "I27", "I28", "I29",
         "I30", "I31", "I32", "I33", "I34", "I35", "I36", "I37", "I38", "I39", "I40",
         "I41", "I42", "I43", "I44", "I45", "I46", "I47", "I48", "I49", "I50", "I51")
hypertension_deaths <- c("I10", "I12", "I15")
cerebrovascular_deaths <- c("I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69")
#add in remaining codes to encompass all I00-I99
complete_list <- c("I14", "I16", "I17", "I18", "I19", "I52", "I53", "I54", "I55", "I56", "I57", "I58", "I59",
                    "I70", "I71", "I72", "I73", "I74", "I75", "I76", "I77", "I78", "I79",
                    "I80", "I81", "I82", "I83", "I84", "I85", "I86", "I87", "I88", "I89",
                    "I90", "I91", "I92", "I93", "I94", "I95", "I96", "I97", "I98", "I99")
all <- c(heart_deaths, hypertension_deaths, cerebrovascular_deaths, complete_list)  #length = 100

#create cvd and alternative info
deaths1$id <- as.character(deaths1$id)
deaths_filt <- deaths1 %>% 
  filter(id %in% pheno$id) #subset to test IDs
  
deaths_filt <- deaths_filt %>% #create new variable for cardiac_death or not
  mutate(cardiac_death = ifelse(
          if_any(everything(), ~. %in% c(all)),
                    "1", "0")) %>% 
    dplyr::select(id,dod_ym, cardiac_death) 

#create column for if death
deaths_filt <- deaths_filt %>%
  mutate(dead = c("1"))

#check numbers
deaths_filt %>% dplyr::filter(dead == "1") %>% dim() 
deaths_filt %>% dplyr::filter(dead == "1" & cardiac_death == "1") %>% dim()
deaths_filt %>% dplyr::filter(dead == "1" & cardiac_death == "0") %>% dim() 
deaths1 %>% dplyr::filter(id %in% pheno$id & if_any(everything(), ~. %in% c(all))) %>% dim() 


#add appt dates
gsappt$id <- as.character(gsappt$id)
pheno_ids <-left_join(pheno, gsappt, by = "id")
dim(pheno_ids)

pheno_select <- pheno_ids %>% 
  dplyr::select(id, age, sex, appt, bmi, pack_years, rank) #%>%   #add in disease Y/N to compare with dates data
  sapply(pheno_select, function(x) sum(is.na(x)))

cvd_merge <- cvd_cases3 %>% dplyr::select(id, dt1_ym, cvd, MI, Isch_stroke, CHD_NOS) #prep cvd df for merging
sapply(cvd_merge, function(x) sum(is.na(x)))  #0 missing values

#add cleaned cvd data to id df
cvd_merge$id <- as.character(cvd_merge$id)
df_cox <- left_join(pheno_select, cvd_merge, by = "id")
sapply(df_cox, function(x) sum(is.na(x)))  #only missing data is for those without cvd 

#add deaths data
deaths_filt$id <- as.character(deaths_filt$id)
df_cox$id <- as.character(df_cox$id)
df_cox <- df_cox %>% left_join(deaths_filt, by = "id")
sapply(df_cox, function(x) sum(is.na(x)))

#check for overlap between cvd cases and deaths & any data entry anomalies
df_cox %>% filter(cvd == "1" & dead == "1") %>% dim()
df_cox %>% filter(cvd == "1" & dead == "1" & cardiac_death == "1") %>% dim() 
df_cox %>% filter(cvd == "1" & dead == "1" & dt1_ym > dod_ym) 
df_cox %>% filter(cvd == "1" & dead == "1" & cardiac_death == "1" & dt1_ym > dod_ym) 

#remove any individuals with prevalent cvd
dim(prevalent_cvd) 
df_cox <- df_cox %>% filter(!id %in% prevalent_cvd$id) 
df_cox %>% dplyr::filter(cvd == "1") %>% dim()
table(df_cox$cvd) 
table(df_cox$cardiac_death)
table(df_cox$dead) 

#set up variables for cox analysis
#create composite outcome, where 1 signifies either cardiac disease or cardiac death and 0 signifies neither
df_cox <- df_cox %>%
  mutate(
    composite_cardiac_outcome = ifelse(
      (coalesce(cvd, "0") == "1" | coalesce(cardiac_death, "0") == "1"), #turns NAs to 0s
      "1",
      "0"
    )
  )
table(df_cox$composite_cardiac_outcome) 


#add in censor date: as this is when event data is from
df_cox <- df_cox %>% mutate(censor_date = c(ym(202310)))
df_cox2 <- df_cox %>% mutate(
    dod_ym = ym(dod_ym) 
)

head(df_cox2)

#check for any data anomalies
# - any CVD events after censor date?
df_cox2 %>% dplyr::filter(dt1_ym > censor_date) 
# - any deaths prior to cvd event
df_cox2 %>% dplyr::filter(dod_ym < dt1_ym) 
df_cox2 %>% dplyr::filter(!is.na(dod_ym)) %>% head()
df_cox2 %>% dplyr::filter(!is.na(dod_ym)) %>% head()
df_cox2 %>% dplyr::filter(dod_ym < censor_date) %>% dim() 
df_cox2 %>% dplyr::filter(cardiac_death == "1" & dod_ym < censor_date) %>% dim() 

df_cox2 <- df_cox2 %>%
  mutate(event_date = pmin(censor_date, dt1_ym, dod_ym, na.rm = TRUE)) #event date is the earliest out of diagnosis/death/censor
sapply(df_cox2, function(x) sum(is.na(x)))

#do some checks
noncvd <- df_cox2 %>% dplyr::filter(dead == "1" & cardiac_death == "0") 
cvd_anddeath <- df_cox2 %>% dplyr::filter(cvd == "1" & cardiac_death == "1") 

#Calculate time to event
 df_cox2$appt2 <- ymd(str_sub(df_cox2$appt, 1,10))

df_cox2 <- df_cox2 %>%
  mutate(tte = difftime(event_date, appt2, units = "days" ))

#check no negative results
df_cox3 <- df_cox2 %>% dplyr::filter(!tte < 0) 
df_cox3 <- df_cox3 %>% dplyr::filter(!id %in% prevalent_cvd$id)
 dim(df_cox3)

df_cox3 <- df_cox3 %>%
  mutate(tte_years = as.numeric(tte)/365.25)

#set status as censored or event

df_cox3 <- df_cox3 %>%
  mutate(status = ifelse(composite_cardiac_outcome == "1" & event_date <= censor_date, "1", "0")) #event on censor date
outcomes <- df_cox2 %>% dplyr::filter(composite_cardiac_outcome == "1") 

  table(df_cox3$status)


#create status columns for IHD, MI, Isch stroke
df_cox3 <- df_cox3 %>%
  mutate(status_IHD = ifelse(CHD_NOS == "1" & event_date <= censor_date, "1", "0"),
        status_MI = ifelse(MI == "1" & event_date <= censor_date, "1", "0"),
        status_stroke = ifelse(Isch_stroke == "1" & event_date <= censor_date, "1", "0")) 

df_cox3 %>% filter(status_IHD == "0" & CHD_NOS == "1") #0
df_cox3 %>% filter(status_MI == "0" & MI == "1") #0
df_cox3 %>% filter(status_stroke == "0" & Isch_stroke == "1")

max(df_cox2$dod_ym, na.rm = TRUE)
max(df_cox2$dt1_ym, na.rm = TRUE)

#___________________________________________________________________________
#2. Filter to bmi >=16 & bmi <=50 

df_cox_filt <- df_cox3 %>%
  dplyr::filter(bmi >= 16 & bmi <= 50 | is.na(bmi)) 
table(df_cox_filt$status) 
table(df_cox_filt$composite_cardiac_outcome) 
table(df_cox_filt$status_MI) 
table(df_cox_filt$status_IHD) 
table(df_cox_filt$status_stroke) 

#Filtered above
# Check for IDs with more than one entry
duplicate_ids <- df_cox_filt %>%
  group_by(id) %>%                  # Group by ID
  filter(n() > 1) %>%               # Filter groups with more than one entry
  distinct(id)
print(duplicate_ids) #nil

#Save out for use in HTN models
fwrite(df_cox_filt, file = "/CVD_basic_coxdf_bmifilt.csv")

#___________________________________________________________________________
#3. Filter to age > 40 and < 69 - to align with SCORE2 range

#check unfiltered cox_df saved
df_cox_filt <- df_cox_filt %>%
  dplyr::filter(age >= 40 & age <= 69) 


#also filter by bmi
df_cox_filt <- df_cox_filt %>%
    dplyr::filter(bmi >= 16 & bmi <= 50 | is.na(bmi))
dim(df_cox_filt)


#count events again
table(df_cox_filt$status)
table(df_cox_filt$composite_cardiac_outcome) 
table(df_cox_filt$status_MI) 
table(df_cox_filt$status_IHD) 
table(df_cox_filt$status_stroke) 


#save out again
fwrite(df_cox_filt, file = "CVVD_basic_coxdf_agebmifilt.csv")
