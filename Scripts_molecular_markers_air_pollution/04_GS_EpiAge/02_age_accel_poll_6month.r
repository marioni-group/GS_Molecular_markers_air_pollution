
library(tidyverse)
library(data.table)
library(broom)

#Load data
ageaccel <- fread("age_acceleration_clocks.csv")
pollution <- fread("6month_avg_exp_foranalysis_scaled.csv")
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
df <- left_join(target2, pheno2, by = "id") 



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

#check for NA values
na_count <-sapply(df2, function(y) sum(length(which(is.na(y)))))
na_count

#Save out file - for demographic information
fwrite(df2, file = "epi_age_demogr_18392_183poll.csv") #This file has bmi filtering only, no other outlier removal


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
fwrite(df3, file = "epi_age_df_sc_filt_183poll.csv")

###Run Models
# Initialize output lists
model_outputs_raw <- list()
tidy_model_outputs <- list()

#set up input variables
clocks <- colnames(ageaccel)
predictors <- clocks[3:9]
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
      paste(pollutant, "~", predictor, "+ age + sex +  (1|id) + log_units + log_bmi + simd2009v2rank + log_smok
                                           + CD4T + CD8T + Bcell + Mono + NK + Eos") 
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
         term = ifelse(term == "DNAmGrimAge", "DNAmGrimAgev2_accel", term))
        

#Calculate confidence intervals
tidy_output2 <- tidy_output2 %>%
  mutate(confint_low = beta - qnorm(0.975)*se_fixed,
         confint_high = beta + qnorm(0.975)*se_fixed)

#Re-order
tidy_output2 <- tidy_output2 %>%
  dplyr::select(pollutant, predictor, model, term, beta, confint_low, confint_high, se_fixed, z, p)

tidy_output3 <- tidy_output2 %>%
  dplyr::filter(predictor %in% c("DNAmGrimAgev2_accel", "PhenoAge_accel", "DunedinPACE_accel"))

#save out full results
fwrite(tidy_output3, file = "lmekin_models_output_183poll.csv")

#save out raw results
saveRDS(model_outputs_raw, file = "lmekin_models_output_raw_183poll.RDS")


#Create summaries
nominal <- 0.05
Pbonf1 <- 0.05/2
Pbonf2 <- 0.05/(2*3)
thresholds <- c(nominal, Pbonf1, Pbonf2)

tidy_output3 %>% dplyr::filter(predictor == term & p < 0.05)
tidy_output4 <- tidy_output3 %>% 
dplyr::filter(predictor == term) %>%
      mutate(significance = ifelse(p < nominal & p >= Pbonf1, "nominal", 
                              ifelse(p < Pbonf1 & p >Pbonf2,  "P < 0.05/2", 
                                ifelse(p < Pbonf2, "P < 0.05/6", "none"))))



install.packages("wesanderson")
library(wesanderson)
colours <- wes_palette("FantasticFox1")
colours <- colours[3:5]

#full model
model4df <- tidy_output4 %>% dplyr::filter(model == "model_4")

#check significant associations at Pbonf1
model4df %>% filter(p < Pbonf1)

#Create plot HRs by predictor 
p2 <- ggplot(data=model4df, aes(y=pollutant, x=beta, xmin=confint_low, xmax=confint_high, colour = significance)) +
  geom_point( size = 2) + 
  geom_errorbarh( height=.1) +
  guides(
    colour = guide_legend("Bonferroni significance")) +
  labs(title='', x='Estimate [95% Confidence Interval]', y = 'Pollutant') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  theme_minimal() + 
  theme(legend.position = "top",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 14)) +
  facet_wrap(~ predictor, ncol = 3)

summary_output <- list()
summary_output <- list()

for (threshold in thresholds) {
  summary_output[[as.character(threshold)]] <- tidy_output3 %>%
    dplyr::filter(p < threshold & predictor == term) %>%
    group_by(predictor, model, pollutant) %>%
    summarise(
      Significant_predictors = paste(unique(term), collapse = ","),
      Pollutant_count = n(),
      
      .groups = "drop"
    )
}

summary_output

fwrite(summary_output, file = "summary_models_2ndgenclocks_183poll.csv")

#summary output by model
summary_output2 <- tidy_output3 %>%
  group_by(model, predictor) %>%
  dplyr::filter(p < Pbonf & term == predictor) %>%
  summarise(
    Significant_pollutants = paste(pollutant, collapse = ","),
    Pollutant_count = n()
  )

fwrite(summary_output2, file = "summary_bymodel_183poll.csv")
