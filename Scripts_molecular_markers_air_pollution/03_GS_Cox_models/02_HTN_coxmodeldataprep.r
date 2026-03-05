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


#Add in HTN tte ***This version has prevalent CVD removed *** 
coxdf <- fread("CVD_basic_coxdf_bmifilt.csv")
colnames(coxdf)

#reload diseases data
diseases <- fread("2025-01-09_diseases.csv")

#filter to IDs in coxdf & to incident cases of hypertension
htn <- diseases %>% 
  filter(incident == "1" & disease == "hypertension" & source == "Secondary_Care")
dim(htn) 
head(htn)

#Select outcome only 
htn2 <- htn %>% 
    mutate(htn = "1",
           htn_dt = lubridate::ym(dt1_ym)) %>%
    dplyr::select(id, htn_dt, htn)

#check for no. prevalent cvd - defined above
head(prevalent_cvd) 
prevalent_cvd %>% filter(id %in% htn2$id) %>% dim() 

#Filter to IDs with pollution
coxdf <- coxdf %>% filter(id %in% annual_df$ID)
dim(coxdf)

coxdf2 <- left_join(coxdf, htn2, by = "id")  
coxdf2 <- coxdf2 %>%
  mutate(status_htn = ifelse(coalesce(htn, "0") == "1", "1", "0"),
          event_htn = pmin(censor_date, htn_dt, dod_ym, na.rm = TRUE),
          tte_htn_days = difftime(event_htn, appt2, units = "days"),
          tte_htn_years = as.numeric(tte_htn_days)/365.25) 

#checks
#check for negative TTEs
coxdf2 %>% filter(status_htn == "1" & event_htn < appt2) %>% dim() #
coxdf3 <- coxdf2

#check for deaths prior to diagnosis
coxdf3 %>% filter(status_htn == "1" & htn_dt > dod_ym) 

table(coxdf3$status_htn) 
table(coxdf3$dead) 
table(coxdf3$dead, coxdf3$status_htn)


#Now create flags for prevalent diabetes & prevalent HTN
flags <- diseases %>% 
  filter(id %in% coxdf2$id) %>%
  mutate(diabetes = ifelse(incident == "0" & disease == "diabetes" & source == "Secondary_Care", "1", "0"),
         hypertension = ifelse(incident == "0" & disease == "hypertension" & source == "Secondary_Care", "1", "0"))

flags2 <- flags %>%
  select(id, diabetes, hypertension) %>%
  filter(id %in% coxdf3$id) %>%
  filter(diabetes == "1" | hypertension == "1")

flags2 <- flags2 %>%
 group_by(id) %>%
 summarise(
  diabetes = max(as.numeric(diabetes)),
  hypertension = max(as.numeric(hypertension))
 )
#merge
coxdf4 <- left_join(coxdf3, flags2, by = "id")  
coxdf4 <- coxdf4 %>%
  mutate(diabetes = coalesce(diabetes, 0),
         hypertension = coalesce(hypertension, 0))  

#check negative ttes and set to prevalent
table(coxdf4$hypertension) 

coxdf4 <- coxdf4 %>% 
  mutate(hypertension = ifelse(tte_htn_years < 0, 1, hypertension))
table(coxdf4$hypertension) 
coxdf4 %>% filter(id %in% prevalent_cvd) #0

#save out ***NB this file still contains those with prevalent hypertension
fwrite(coxdf4, file = "CVD_htn_basic_coxdf_bmifilt.csv")  
