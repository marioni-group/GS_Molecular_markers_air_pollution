#Pollution - Dementia 

#. 1. Set-up incident diseases of interest - generate yes(1) or no(0) column (Dementia: diagnosis or diagnosis at death)
#. 2. Remove prevalent disease (individuals diagnosed before onset of study)
#. 3. Filter to age at diagnosis >65 or age at censoring >65 (to remove anomalies e.g. early onset dementia)


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

#Load Data
data_path <- here::here("")
annual_df <- fread("annual_avg_exp_foranalysis_ntr.csv") 
dementia_cases <- fread("dementia_case_update_20241024.csv") #Disease phenotypes. Censor date will be October '24
deaths1 <- fread(here("2024-10-24_deaths.csv")) #CVD deaths
gsappt <- fread(here("GS_appt.txt")) #appointment dates
pheno <- readRDS("GS_phenos_internal_23Aug2023_REM.rds") #other covariates
simd <- fread("2024-10-18_simd09_all_categories.csv") 
simd <- simd %>% dplyr::rename(ID = id)

#________________________________________________________________________________________
#1. Generate outcome df for individuals in test set: dementia cases or not

#Filter pheno to ids with pollution data
pheno$id <- as.character(pheno$id)
pheno <- pheno %>% dplyr::filter(id %in% annual_df$ID)

#subset only to columns of interest for creating the cox df
pheno_cox <- pheno %>% dplyr::select(id, age, sex, dob_ym, y_appt)

#subset disease phenotypes to ids in pheno
dementia_cases <- dementia_cases %>% dplyr::filter(.,id %in% pheno$id)
length(unique(dementia_cases$id))

#_________________________________________________________________________________________
#2. Remove prevalent cases
#Look at prevalent disease
class(dementia_cases$appt)
prevalent_dementia <- dementia_cases %>%
  dplyr::filter(dt1_dementia < appt) 
length(unique(prevalent_dementia$id))

#Filter out prevalent cases
dementia_cases2 <- dementia_cases %>%
  dplyr::filter(!id %in% prevalent_dementia$id)
dim(dementia_cases2)
length(unique(dementia_cases2$id))

#___________________________________________________________________________
#3. Filter to age at diagnosis > 65 or age at censoring >65

#Check for any cases in individuals < 65
dementia_cases2$id <- as.character(dementia_cases2$id)
dementia_cases3 <- left_join(dementia_cases2, pheno_cox, by = "id")
# dim(dementia_cases3)

#create column for age at diagnosis
dementia_cases4 <- dementia_cases3 %>% 
  mutate(
    dt1_dementia = ym(dt1_dementia),
    appt = ym(appt))

dementia_cases4 <- dementia_cases4 %>%
  mutate(
    age_diagnosis = age + as.numeric(difftime(dt1_dementia,appt, units = "weeks"))/52.1775)

agefilt_d <- dementia_cases4 %>% dplyr::filter(age_diagnosis >= 65)
length(unique(agefilt_d$id)) 
agefilt_d %>% dplyr::filter(dt1_dementia == appt) 
agefilt_d %>% dplyr::filter(dt1_dementia < appt) 
agefilt_d %>% dplyr::filter(dt1_dementia > appt) %>% dim() 

#retain earliest diagnosis date where there are multiple entries
agefilt_unique <- agefilt_d %>%
  arrange(dt1_dementia) 

agefilt_unique <- agefilt_d %>%
  dplyr::group_by(id) %>% #for each id, filter to earliest diagnosis date
  mutate(is_death = any(source == "deaths" & dt1_dementia == min(dt1_dementia))) %>%
   filter(dt1_dementia == min(dt1_dementia)) %>%
   ungroup() %>%
   distinct(., id, .keep_all = TRUE)
head(agefilt_unique)
dim(agefilt_unique) 
#This df contains data for individuals with a dementia diagnosis

table(agefilt_unique$is_death)


#clean up agefilt_unique
agefilt_unique <- agefilt_unique %>% 
 dplyr::select(-c(age, sex, dob_ym, y_appt, ad, vd, 
                    dlb, ftd, n, dt1_ad, dt1_vd, dt1_dlb, 
                    dt1_ftd, dt1_n))

#_____________________________________________
#Now create cox df for all

##Prep input variables

#create cox dataframe with diagnosis, pheno and pollution data by ID
pheno_cox$id <- as.character(pheno_cox$id) #Pheno_cox contains variables relevant to creating cox outcomes
 
cox_df <- dplyr::left_join(pheno_cox, agefilt_unique, by = "id", keep = FALSE)
colnames(cox_df)

#add in censor date
cox_df2 <- cox_df %>%
  mutate(
  censor_date = ym(c("2024-10")))

#add in age at censor date
#First need to add in appt for everyone
head(gsappt)
gsappt2 <- gsappt %>%
  mutate(id = as.character(id),
         appt2 = ymd(str_sub(appt, 1,10)))

cox_df3 <- left_join(cox_df2, gsappt2, by = "id")

cox_df3 <- cox_df3 %>%
  mutate(
    age_censor = age + as.numeric(difftime(censor_date,appt2, units = "weeks"))/52.1775)  

dim(cox_df3) 
cox_df4 <- cox_df3 %>% dplyr::filter(age_censor >= 65)
dim(cox_df4) 
head(cox_df4)
# # #check for any missing values
sapply(cox_df4, function(x) sum(is.na(x))) 

#clean up df
cox_df4 <- cox_df4 %>% 
  dplyr::select(-c(appt.x, appt.y))

#Add in deaths
head(deaths1)
deaths1$id <- as.character(deaths1$id)
cox_df5 <- dplyr::left_join(cox_df4, deaths1, by = "id")
dim(cox_df5)


#Make sure any prevalent dementia cases are removed
cox_df5 <- cox_df5 %>% dplyr::filter(!id %in% prevalent_dementia$id)
dim(cox_df5) 


#Check for any data entry anomalies
# - any dementia events after censor date?
cox_df5 %>% dplyr::filter(dt1_dementia > censor_date) 
# - any dementia cases before appt date
cox_df5 %>% dplyr::filter(dt1_dementia < appt2) 


#___________________________________________________
#Add in outcome variables to be used in cox df
#1. event date (dementia_diagnosis/dementia_related_death/non-dementia death/censorship)
#2. TTE (in years)
#3. Status

#.1 Event date
#check formats of relevent dates
str(cox_df5)
cox_df5$dod_ym <- ym(cox_df5$dod_ym)

cox_df5 <- cox_df5 %>%
  mutate(event_date = pmin(censor_date, dt1_dementia, dod_ym, na.rm = TRUE)) #event date is the earliest out of diagnosis/death/censor
sapply(cox_df5, function(x) sum(is.na(x)))

#.2 TTE
cox_df5 <- cox_df5 %>%
  mutate(tte = difftime(event_date, appt2, units = "days" ))
#check no negative results
cox_df5 %>% dplyr::filter(tte < 0) 
#create time to event in years
cox_df5 <- cox_df5 %>%
  mutate(tte_years = as.numeric(tte)/365.25)

#.3 Status
#set status as censored or event
cox_df5$dementia <- replace(cox_df5$dementia, is.na(cox_df5$dementia), "0")
cox_df5 <- cox_df5 %>%
  mutate(status = ifelse(dementia == "1" & dt1_dementia == event_date & event_date < censor_date, "1", 
                     ifelse(dementia == "1" & is_death == "TRUE" & event_date < censor_date, "1", "0")))
cox_df5$status <- as.numeric(cox_df5$status) 

cox_df5 %>% dplyr::filter(dementia == "1") %>% dim() #323
table(cox_df5$status)

#check for missing values
sum(is.na(cox_df5$tte_years))
summary(cox_df5$tte_years)

#save out
fwrite(cox_df5, file = "dementia_coxdf_10761.csv")
