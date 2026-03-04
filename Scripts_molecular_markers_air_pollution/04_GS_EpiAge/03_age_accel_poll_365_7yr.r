#From: /Cluster_Filespace/Marioni_Group/Josie/Proteins/scripts/pollution_analysis/EpiAge/april_2025/age_accel.poll.r

library(tidyverse)
library(data.table)
library(broom)

#Load data
# clock_data<- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/gs_clocks.csv")
ageaccel <- fread("age_acceleration_clocks.csv")
pollution <- fread("annual_avg_exp_foranalysis_scaled.csv")
target <- readRDS("GS20k_Targets_18869.rds")
pheno <- readRDS("GS_phenos_internal_23Aug2023_REM.rds")
simd <- fread("2024-10-18_simd09_all_categories.csv")

#Make a kinship matrix
library(kinship2)
#read in pedigree file
ped1 <- read.csv("2023-03-20_pedigree.csv")
#input missing fam IDs
ped2 = data.frame(famid=c(4091,4384), volid=c(103027, 144865), father=c(0,0), mother=c(0,0), sex= c("F", "F"))

#bind together
ped = rbind(ped1, ped2)
#
#make kinship matrix
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin)

#Items to iterate over
clocks <- colnames(ageaccel)
predictors <- clocks[3:9]
pollutants <- colnames(pollution)[2:9]

#create covariate df
pheno2 <- pheno %>% dplyr::select(id, age, sex, bmi, units, pack_years)
simdrank <- simd %>% dplyr::select(id, simd2009v2rank)
pheno2 <- left_join(pheno2, simdrank, by = "id")
target2 <- target %>% dplyr::select(Sample_Name, Sample_Sentrix_ID, Batch, Bcell, CD4T, CD8T, Eos, Mono, NK)
target2 <- target2 %>% dplyr::rename(id = Sample_Name)
pheno2$id <- as.character(pheno2$id)
df <- left_join(target2, pheno2, by = "id") #18869


#add pollution data & clock data
pollution$id <- as.character(pollution$ID)

#filter df to ids in pollution data
df <- df %>% dplyr::filter(id %in% pollution$id) 
meth_ids <- fread("2024-10-31_meth_ids.csv")
df_check <- df %>% dplyr::filter(id %in% meth_ids$V1) 

#filter ageaccel to pollution ids
ageaccel <- ageaccel %>% filter(Sample_Sentrix_ID %in% df$Sample_Sentrix_ID) 

#add 
df1 <- left_join(df, ageaccel, by = "Sample_Sentrix_ID") 
pollution <- pollution %>% filter(id %in% df1$id) 
df2 <- left_join(df1, pollution, by = "id") 


#filter by BMI
dim(df2) 
df2 <- df2 %>% dplyr::filter(bmi >=16 & bmi <= 50)
dim(df2) 

#tidy up dataframe
df2 <- df2 %>% select(-age.y) %>%
  rename(age = age.x) 
df2 <- df2 %>% select(-c(ID))
df2 <- df2 %>% select(c(id, Sample_Sentrix_ID, Batch, sex, everything())) 
df2 <- df2 %>% dplyr::rename(DNAmGrimAgev2_accel = "DNAmGrimAge.1_accel") #Rename to avoid confusion
fwrite(df2, file = "demographics_raw_epiage.csv")

#Summarise demographics for supplementary
#Create summary table
#Extract demographic info from non-imputed covariates
library(arsenal)
table <- tableby(sex ~ age + bmi + pack_years + units + simd2009v2rank, 
                  data = df2, 
                  test = FALSE, 
                  total = TRUE,
                  control = tableby.control(
                  numeric.stats = c("mean", "sd", "median", "q1q3", "Nmiss2")))
tab1 <- as.data.frame(summary(table))
fwrite(tab1, file = "summary_demo_epiage_18392.csv")

#check for NA values
na_count <-sapply(df2, function(y) sum(length(which(is.na(y)))))
na_count


#Save out file - for demographic information
fwrite(df2, file = "epi_age_demogr_18392.csv")


#Scale numeric variables, except pollutants (already scaled)
#First transform skewed variables: units, pack_years, bmi
str(df2)
df2$units <- as.numeric(df2$units)
df2$pack_years[df2$pack_years == "NA"] <- NA 
df2$pack_years <- as.numeric(df2$pack_years)
df2$log_units <- log10(df2$units + 1)
df2$log_smok <- log10(df2$pack_years + 1)
df2$log_bmi <- log10(df2$bmi)

#re-order to simplify scaling
df2 <- df2 %>% dplyr::select(id, Sample_Sentrix_ID, Batch, PM25, NO2, NO, O3, PM10, SO2, NO3_C, NO3_F, everything())
df2[12:33] <- df2[12:33] %>% purrr::modify_if(is.numeric, scale)
sapply(df2, function(y) sum(length(which(is.na(y)))))

#Check if any rows/columns have missing data > 0.4
max_missing <- 0.4
rows_with_high_missing <- which(rowMeans(is.na(df2)) > max_missing) #0
cols_with_high_missing <- which(colMeans(is.na(df2)) > max_missing) #0

#Impute missing data
library(VIM)
set.seed(12)
df3 <- kNN(df2, k = 5)

#Save out prepped df for models
fwrite(df3, file = "epi_age_df_sc_filt.csv")
df3 <- fread("epi_age_df_sc_filt.csv")

#remove imputed variables
df3 <- df3 %>% select(!contains("imp")) #remove columns added by knn step

###Run Models
# Initialize output lists
model_outputs_raw <- list()
tidy_model_outputs <- list()

#set up input variables
clocks <- colnames(ageaccel)
predictors <- clocks[6:9] #select 2nd gen clocks only
predictors <- predictors[-3] #remove GrimAge version 1
pollutants <- colnames(pollution)[2:9]

# Load function to Extract Lmekin Results
extract_lmekin_table <- function (mod){
  beta <- fixef(mod) #$fixed is not needed
  nvar <- length(beta)
  frail <- ranef(mod)
  nfrail <- length(frail)
  variance <- vcov(mod)
  varr_random <- VarCorr(mod)
  se_fixed <- sqrt(diag(variance))
  z<- round(beta/se_fixed, 2)
  p<- 1 - pchisq((beta/se_fixed)^2, 1)
  table=data.frame(cbind(beta,se_fixed, z,p))
  return(table)
}

#Load coxme package
library(coxme)


# Iterate over pollutants and predictors
for (pollutant in pollutants) {
  for (predictor in predictors) {
    print(pollutant)
    print(predictor)
    #model list
    formula_list <- list(
      paste(pollutant, "~", predictor, "+ age + sex + (1|id)"),
      paste(pollutant, "~", predictor, "+ age + sex  + (1|id) + log_units + log_bmi + simd2009v2rank"), 
      paste(pollutant, "~", predictor, "+age + sex + (1|id) + log_units + log_bmi + simd2009v2rank + log_smok"),
      paste(pollutant, "~", predictor, "+ age + sex +  (1|id) + log_units + log_bmi + simd2009v2rank + 
                                        log_smok + CD4T + CD8T + Bcell + Mono + NK + Eos")
    )
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

for (formula in formula_list) {
print(formula)
model <- lmekin(formula, 
  varlist = kin_model*2,
   data = df3, 
   na.action = na.exclude) 
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, predictor, model_index, sep = "_")
model_outputs_raw[[model_output_name]] <- model
tidy_model_outputs[[model_output_name]] <- extract_lmekin_table(model) %>% mutate(model = model_name,
                                pollutant = pollutant,
                                predictor = predictor)
   # Increment model index
      model_index <- model_index + 1
}
tidy_output <- bind_rows(tidy_model_outputs)

  }
}
head(tidy_output)
head(tidy_output$model)
dim(tidy_output)

#Tidy up output
tidy_output2 <- tidy_output %>%
  dplyr::mutate(term = rownames(tidy_output))
tidy_output2$term <- gsub("\\...*", "", tidy_output2$term)
#rename DNAmGrmAge.1_accel so it doesn't cause problems
tidy_output2 <- tidy_output2 %>%
  mutate(predictor = ifelse(predictor == "DNAmGrimAge.1_accel", "DNAmGrimAgev2_accel", predictor),
         term = ifelse(term == "DNAmGrimAge", "DNAmGrimAgev2_accel", term)) #because the original input for DNAmGrimAge had a .1 in it, the term got renamed incorrectly
        

#Calculate confidence intervals
tidy_output2 <- tidy_output2 %>%
  mutate(confint_low = beta - qnorm(0.975)*se_fixed,
         confint_high = beta + qnorm(0.975)*se_fixed)

#Re-order
tidy_output2 <- tidy_output2 %>%
  dplyr::select(pollutant, predictor, model, term, beta, confint_low, confint_high, se_fixed, z, p)


#save out full results
fwrite(tidy_output2, file = "lmekin_models_2ndgen_output.csv")

#save out raw results
saveRDS(model_outputs_raw, file = "lmekin_models_2ndgen_output_raw.RDS")


#________________________________________________________
#Re-run with 7-year pollution data
#__________________________________________________________


poll7yr <- fread("sevenyr_foranalysis_sc.csv")
poll7yr <- poll7yr %>% rename(id = ID)
dim(df3)
df4 <- df3 %>% select(!all_of(pollutants)) #remove 365-day data
colnames(df4)
df4 <- left_join(df4, poll7yr, by = "id") #now contains 7-year data

#run models as above

for (pollutant in pollutants) {
  for (predictor in predictors) {
    print(pollutant)
    print(predictor)
    #model list
    formula_list <- list(
      paste(pollutant, "~", predictor, "+ age + sex + (1|id)"),
      paste(pollutant, "~", predictor, "+ age + sex  + (1|id) + log_units + log_bmi + simd2009v2rank"), 
      paste(pollutant, "~", predictor, "+age + sex + (1|id) + log_units + log_bmi + simd2009v2rank + log_smok"),
      paste(pollutant, "~", predictor, "+ age + sex +  (1|id) + log_units + log_bmi + simd2009v2rank + log_smok +
                                           CD4T + CD8T + Bcell + Mono + NK + Eos")
    )
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

for (formula in formula_list) {
print(formula)
model <- lmekin(formula, 
  varlist = kin_model*2,
   data = df4, 
   na.action = na.exclude) 
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, predictor, model_index, sep = "_")
model_outputs_raw[[model_output_name]] <- model
tidy_model_outputs[[model_output_name]] <- extract_lmekin_table(model) %>% mutate(model = model_name,
                                pollutant = pollutant,
                                predictor = predictor)
   # Increment model index
      model_index <- model_index + 1
}
tidy_output <- bind_rows(tidy_model_outputs)

  }
}

#Tidy up output
tidy_7yr <- tidy_output %>%
  dplyr::mutate(term = rownames(tidy_output))
tidy_7yr$term <- gsub("\\...*", "", tidy_7yr$term)
#rename DNAmGrmAge.1_accel so it doesn't cause problems
tidy_7yr <- tidy_7yr %>%
  mutate(predictor = ifelse(predictor == "DNAmGrimAge.1_accel", "DNAmGrimAgev2_accel", predictor),
         term = ifelse(term == "DNAmGrimAge", "DNAmGrimAgev2_accel", term)) #because the original input for DNAmGrimAge had a .1 in it, the term got renamed incorrectly
        

#Calculate confidence intervals
tidy_7yr <- tidy_7yr%>%
  mutate(confint_low = beta - qnorm(0.975)*se_fixed,
         confint_high = beta + qnorm(0.975)*se_fixed)

#Re-order
tidy_7yr <- tidy_7yr %>%
  dplyr::select(pollutant, predictor, model, term, beta, confint_low, confint_high, se_fixed, z, p)

fwrite(tidy_7yr, file = "lmekin_models_2ndgen_7yrpoll.csv")
tidy_7yr <- fread("lmekin_models_2ndgen_7yrpoll.csv")

#prepare for joining/plotting
poll7yr <- tidy_7yr %>%
  filter(predictor == "PhenoAge_accel" | predictor == "DunedinPACE_accel" | predictor == "DNAmGrimAgev2_accel") %>%
  mutate(significance = ifelse(p < 0.05 & p > 0.05/2, "*", 
                          ifelse(p < 0.05/2 & p > 0.05/6, "**", 
                            ifelse( p < 0.05/6, "***", "")))) %>%
  dplyr::filter(term == predictor) %>%
  mutate(data_type = "7_yr") %>%
  filter(model == "model_1" | model == "model_4")


#_______________________________________________________________________________________
#plot for publ - keep 2nd gen clocks 


secondgen <- fread("filt_models_output.csv")


#Plot
#colour by model, label by significance threshold
thresholds <- c(0.05, 0.05/2, 0.05/6) #thresholds /2 for pollutant PCs, /6 for pollutant PCs * clocks
secondgen <- secondgen %>% 
  mutate(significance = ifelse(p < 0.05 & p > 0.05/2, "*", 
                          ifelse(p < 0.05/2 & p > 0.05/6, "**", 
                            ifelse( p < 0.05/6, "***", ""))
      
  ))
secondgen <- secondgen %>% mutate(data_type = "365_day")

plotdf <- secondgen %>% dplyr::filter(term == predictor)

basic_plot <- ggplot(data=plotdf, aes(y=pollutant, x=beta, xmin=confint_low, xmax=confint_high, colour = model)) +
 geom_point(position = position_dodge(width = 0.8), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.8), height=.1) +
  geom_text(aes(label=significance), position=position_dodge2(width=0.8), vjust = -0.75, size= 2, colour = "black") +
  guides(
    colour = guide_legend("Model")) +
  labs(title='365-day pollution exopsure and biological age acceleration', x='Estimate [95% Confidence Interval]', y = 'Pollutant') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(-0.05, 0.1)) +
  theme_minimal() + 
  theme(legend.position = "top",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 14)) +
  facet_wrap(~ predictor, ncol = 5)
ggsave(basic_plot, file = "clock_poll_2ndgen_allmodels.pdf", height = 30, width = 30, units = "cm")

#_____________________________________________
#add in 7yr, 6month data and 3 month data
#_____________________________________________


#6 month data
poll183 <- fread("lmekin_models_output_183poll.csv")
poll183 <- poll183 %>%
  filter(predictor == "PhenoAge_accel" | predictor == "DunedinPACE_accel" | predictor == "DNAmGrimAgev2_accel") %>%
  mutate(significance = ifelse(p < 0.05 & p > 0.05/2, "*", 
                          ifelse(p < 0.05/2 & p > 0.05/6, "**", 
                            ifelse( p < 0.05/6, "***", "")))) %>%
  dplyr::filter(term == predictor) %>%
  mutate(data_type = "183_day") %>%
  filter(model == "model_1" | model == "model_4")

plot183 <- ggplot(data=poll183, aes(y=pollutant, x=beta, xmin=confint_low, xmax=confint_high, colour = model)) +
 geom_point(position = position_dodge(width = 0.75), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.1) +
  geom_text(aes(label=significance), position=position_dodge2(width=0.75), vjust = -0.75, size= 2, colour = "black") +
  guides(
    colour = guide_legend("Model")) +
  labs(title='183-day pollution exposure and biological age acceleration', x='Estimate [95% Confidence Interval]', y = 'Pollutant') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  theme_minimal() + 
  theme(legend.position = "top",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 14)) +
  facet_wrap(~ predictor, ncol = 5)
ggsave(basic_plot, file = "clock_poll_2ndgen_183day.pdf", height = 20, width = 30, units = "cm")

poll183 <- poll183 %>% select(-data)
poll183 <- poll183 %>% rename(se = se_fixed)
poll183 <- poll183 %>% rename(statistic = z)
secondgen <- secondgen %>% select(-data)
poll7yr <- poll7yr %>%
  rename(se = se_fixed,
        statistic = z)
poll7yr <- poll7yr %>% filter(model == "model_4")


plot183 <- ggplot(data=poll183, aes(y=pollutant, x=beta, xmin=confint_low, xmax=confint_high, colour = model)) +
 geom_point(position = position_dodge(width = 0.75), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.1) +
  geom_text(aes(label=significance), position=position_dodge2(width=0.75), vjust = -0.75, size= 2, colour = "black") +
  guides(
    colour = guide_legend("Model")) +
  labs(title='183-day pollution exposure and biological age acceleration', x='Estimate [95% Confidence Interval]', y = 'Pollutant') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  theme_minimal() + 
  theme(legend.position = "top",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 14)) +
  facet_wrap(~ predictor, ncol = 5)
ggsave(basic_plot, file = "clock_poll_2ndgen_183day.pdf", height = 20, width = 30, units = "cm")

poll183 <- poll183 %>% select(-data)
poll183 <- poll183 %>% rename(se = se_fixed)
poll183 <- poll183 %>% rename(statistic = z)
secondgen <- secondgen %>% select(-data)
poll7yr <- poll7yr %>%
  rename(se = se_fixed,
        statistic = z)
poll3mnth <- poll3mnth %>% 
  rename(se = se_fixed,
        statistic = z)

dfc <- rbind(poll183, secondgen, poll7yr)
dfc <- dfc %>% filter(model == "model_4") %>% filter(term == predictor)
fwrite(dfc, file = "model4_diffduratpoll_results.csv")

#summarise dfc
summary <- dfc %>%
 filter(!significance == "", !significance == "*") %>%
 group_by(data_type, pollutant, significance) %>%
 summarise(
  duration = data_type,
  pollutant = pollutant,
  significance = significance,
  significant_clocks = paste(predictor)
 )
summary2 <- dfc %>%
  filter(!significance %in% c("", "*")) %>%
  group_by(data_type) %>%
  summarise(
    pollution_output = paste(unique(pollutant), collapse = ", "),
    clocks = paste(unique(predictor), collapse = ", "),
    .groups = 'drop'
  )
fwrite(summary2, file = "summary_model4_diffduratpoll_results.csv")

summary_output3 <- dfc %>%
filter(!significance %in% c("", "*")) %>%
  group_by(predictor, data_type) %>%
  summarise(
    pollutant = paste(unique(pollutant), collapse = ", "),
    .groups = "drop"
  )
fwrite(summary_output3, file = "summary_model4bypredictor_diffduratpoll_results.csv")

combplot <- ggplot(data=dfc, aes(y=pollutant, x=beta, xmin=confint_low, xmax=confint_high, colour = data_type)) +
 geom_point(position = position_dodge(width = 0.8), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.8), height=.1) +
  geom_text(aes(label=significance), position=position_dodge2(width=0.8), vjust = -0.55, size= 2, colour = "black") +
  guides(
    colour = guide_legend("Model")) +
  labs(title='Pollution exposure of different durations and biological age acceleration', x='Estimate [95% Confidence Interval]', y = 'Pollutant') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  theme_minimal() + 
  theme(legend.position = "top",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 14)) +
  facet_wrap(~ predictor, ncol = 5)
ggsave(combplot, file =  "clock_poll_2ndgen_durationcomp.pdf", height = 20, width = 30, units = "cm")
