#Look at associations of pollution PCs & Health outcomes
#Start with CVD-related outcomes
#need to be in micromamba environment to use pcatools

#use micromamba environment: in terminal run:
eval "$(micromamba shell hook --shell bash)"
micromamba activate R
R


#load libraries
library(tidyverse)
library(data.table)
library(coxme)
library(survival)
library(survminer)
library(knitr)
library(missMethods)
library(car)
library(caret)
library(pscl)
library(PCAtools)

#Load PCs
pollpca <- readRDS(file = "/pollution_pca.RDS")
rotations <- pollpca$rotated
dim(rotations)

#keep first 2 PCs
pc <- rotations %>% select(PC1, PC2) %>% mutate(id = rownames(.))

#load coxdf 
cox_df <- fread("/cvd_cox_df_agebmifilt_sc_imp_070725.csv")

#load kinship matrix
kin_model <- readRDS("/Dementia/kin_model_0901.rds")

#function to extract results:
#function to extract coxme results (tidy doesn't work on coxme format)
extract_coxme_results <- function(model, pollutant, model_name) {
  fixed_effects <- summary(model)$coefficients   #fixed effects
  conf_intervals <- exp(confint(model)) #confidence intercals
  num_terms <- nrow(fixed_effects)
  results_df <- data.frame( 
    pollutant = rep(pollutant, num_terms),        #results df
    term = rownames(fixed_effects),
    estimate = exp(fixed_effects[, "coef"]),
    std.error = fixed_effects[, "se(coef)"],
    statistic = fixed_effects[,"z"],
    p.value = fixed_effects[,"p"],
    conf.low = conf_intervals[, 1],
    conf.high = conf_intervals[1:num_terms, 2],
    row.names = NULL
  ) 
  return(results_df)
}

#function to get AIC-like metric
extract_coxme_aic <- function(model, pollutant, model_name) {
  logLik_value <- logLik(model)
  k <- attr(logLik_value, "df")
  aic_value <- -2 * as.numeric(logLik_value) + 2 * k
  aic_results <- data.frame(
    pollutant = pollutant,
    model = model_name,
    aic_value = aic_value
  )
  return(aic_results)
}

#Ischaemic stroke
# #Set up variables to iterate over
pollutants1 <- c("PC1", "PC2")#pollutants stratified into high vs. low


#join PCs to coxdf
cox_df$id <- as.character(cox_df$id)
coxdf2 <- left_join(cox_df, pc, by = "id")


formula <- as.formula("Surv(tte_years, status_stroke) ~ PC2 + age + as.factor(sex)  + 
                      (1|id) + logunits + logbmi + logsmok + simd2009v2rank + as.factor(diabetes)
                       + Total_cholesterol + HDL_cholesterol + as.factor(hypertension)")
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
summary(model)


#Add in categorised pollutants
for (pollutant in pollutants1) {
  coxdf2 <- coxdf2 %>%
    mutate(
      !!paste0(pollutant, "b") := cut(
        .data[[pollutant]], 
        breaks = quantile(.data[[pollutant]], probs = 0:2 / 2, na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2")
      ))}

formula <- as.formula("Surv(tte_years, status_stroke) ~ PC2b + age + as.factor(sex) +
                           (1|id) + logunits + logbmi + logsmok + simd2009v2rank + 
                           as.factor(diabetes) + Total_cholesterol + HDL_cholesterol + 
                           as.factor(hypertension)")
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
summary(model)


#make sure outcomes have NAs as 0 
coxdf2 <- coxdf2 %>%
  mutate(status_stroke = coalesce(status_stroke, 0),
         status_IHD = coalesce(status_IHD, 0),
         status_MI = coalesce(status_MI, 0))

assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors
pollutants2 <- c("PC1b", "PC2b")

for (pollutant in pollutants2) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status_stroke) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status_stroke) ~", pollutant, "+ age + as.factor(sex)  + (1|id) +
                                                              logunits + logbmi + logsmok + simd2009v2rank"),
      paste("Surv(tte_years, status_stroke) ~", pollutant, "+ age + as.factor(sex)  + (1|id) +
                                                                 logunits + logbmi + logsmok + simd2009v2rank 
                                                                 + as.factor(diabetes) + Total_cholesterol + 
                                                                 HDL_cholesterol + as.factor(hypertension)")
    )

    
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

   #initialise lists for models
   model_results <- list()
   assump_results <- list()
   aic_list <- list()
for (formula in formula_list) {
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>%
                                               mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% 
                                as.data.frame() %>% 
                              mutate(model = model_name, pollutant = pollutant)
aic_list[[model_name]] <- extract_coxme_aic(model, pollutant, model_name)
   # Increment model index
      model_index <- model_index + 1     
}
assump_test[[pollutant]] <- bind_rows(assump_results) 
cox_results[[pollutant]] <- bind_rows(model_results) 
aic_results[[pollutant]] <- bind_rows(aic_list)
  }

cox_results_df <- bind_rows(cox_results)
assump_test_df <- bind_rows(assump_test)
assump_test_df$term <- rownames(assump_test_df)
aic_results_df <- bind_rows(aic_results)

fwrite(cox_results_df, file = "/Pollution_PCs/coxme_agebmipollcat_stroke_0707.csv")
fwrite(assump_test_df, file = "/Pollution_PCs/cox_assump_agebmipollcat__stroke_0707.csv")
fwrite(aic_results_df, file = "/Pollution_PCs/cox_aic_agebmipollcat_stroke_0707.csv")


#IHD
assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants2) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status_IHD) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status_IHD) ~", pollutant, "+ age + as.factor(sex)  + (1|id)
                                                             + logunits + logbmi + logsmok + 
                                                                  simd2009v2rank"),
      paste("Surv(tte_years, status_IHD) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                          logunits + logbmi + logsmok + 
                                                          simd2009v2rank + as.factor(diabetes) + 
                                                          Total_cholesterol + HDL_cholesterol + 
                                                          as.factor(hypertension)")
    )

    
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

   #initialise lists for models
   model_results <- list()
   assump_results <- list()
   aic_list <- list()
for (formula in formula_list) {
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>%
                                                     mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% as.data.frame() %>%
                                                         mutate(model = model_name, pollutant = pollutant)
aic_list[[model_name]] <- extract_coxme_aic(model, pollutant, model_name)
   # Increment model index
      model_index <- model_index + 1     
}
assump_test[[pollutant]] <- bind_rows(assump_results) 
cox_results[[pollutant]] <- bind_rows(model_results) 
aic_results[[pollutant]] <- bind_rows(aic_list)
  }

cox_results_df <- bind_rows(cox_results)
assump_test_df <- bind_rows(assump_test)
assump_test_df$term <- rownames(assump_test_df)
aic_results_df <- bind_rows(aic_results)

fwrite(cox_results_df, file = "/Pollution_PCs/coxme_agebmipollcat_ihd_0707.csv")
fwrite(assump_test_df, file = "/Pollution_PCs/cox_assump_agebmipollcat_ihd_0707.csv")
fwrite(aic_results_df, file = "/Pollution_PCs/cox_aic_agebmipollcat_ihd_0707.csv")


# #MI
assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants2) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status_MI) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status_MI) ~", pollutant, "+ age + as.factor(sex)  + (1|id) 
                                                        + logunits + logbmi + logsmok + 
                                                        simd2009v2rank"),
      paste("Surv(tte_years, status_MI) ~", pollutant, "+ age + as.factor(sex)  + (1|id)
                                                             + logunits + logbmi + logsmok + 
                                                             simd2009v2rank + as.factor(diabetes) + 
                                                             Total_cholesterol + HDL_cholesterol + 
                                                             as.factor(hypertension)")
    )

    
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

   #initialise lists for models
   model_results <- list()
   assump_results <- list()
   aic_list <- list()
for (formula in formula_list) {
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>% mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% as.data.frame() %>% mutate(model = model_name, pollutant = pollutant)
aic_list[[model_name]] <- extract_coxme_aic(model, pollutant, model_name)
   # Increment model index
      model_index <- model_index + 1     
}
assump_test[[pollutant]] <- bind_rows(assump_results) 
cox_results[[pollutant]] <- bind_rows(model_results) 
aic_results[[pollutant]] <- bind_rows(aic_list)
  }

cox_results_df <- bind_rows(cox_results)
assump_test_df <- bind_rows(assump_test)
assump_test_df$term <- rownames(assump_test_df)
aic_results_df <- bind_rows(aic_results)

fwrite(cox_results_df, file = "/Pollution_PCs/coxme_agebmipollcat_mi_0707.csv")
fwrite(assump_test_df, file = "/Pollution_PCs/cox_assump_agebmipollcat_mi_0707.csv")
fwrite(aic_results_df, file = "/Pollution_PCs/cox_aic_agebmipollcat_mi_0707.csv")

# #all cardiac
assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants2) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex)  + (1|id) +
                                                         logunits + logbmi + logsmok + simd2009v2rank"),
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex)  + (1|id) +
                                                               logunits + logbmi + logsmok + 
                                                               simd2009v2rank + as.factor(diabetes) +
                                                                Total_cholesterol + HDL_cholesterol +
                                                                 as.factor(hypertension)")
    )

    
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

   #initialise lists for models
   model_results <- list()
   assump_results <- list()
   aic_list <- list()
for (formula in formula_list) {
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>%
                                                               mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% 
                                    as.data.frame() %>% 
                                mutate(model = model_name, pollutant = pollutant)
aic_list[[model_name]] <- extract_coxme_aic(model, pollutant, model_name)
   # Increment model index
      model_index <- model_index + 1     
}
assump_test[[pollutant]] <- bind_rows(assump_results) 
cox_results[[pollutant]] <- bind_rows(model_results) 
aic_results[[pollutant]] <- bind_rows(aic_list)
  }

cox_results_df <- bind_rows(cox_results)
assump_test_df <- bind_rows(assump_test)
assump_test_df$term <- rownames(assump_test_df)
aic_results_df <- bind_rows(aic_results)

fwrite(cox_results_df, file = "/Pollution_PCs/coxme_agebmipollcat_allcardiac.csv")
fwrite(assump_test_df, file = "/Pollution/Cox_models/Pollution_PCs/cox_assump_agebmipollcat_allcardiac.csv")
fwrite(aic_results_df, file = "/Pollution/Cox_models/Pollution_PCs/cox_aic_agebmipollcat_allcardiac.csv")


#_______________new coxdf for htn analysis_____________#

cox_df <- fread("htn_cox_df_bmifilt_sc_imp_070725.csv")
cox_df$id <- as.character(cox_df$id)

coxdf2 <- left_join(cox_df, pc, by = "id")


coxdf2 <- coxdf2
for (pollutant in pollutants1) {
  coxdf2 <- coxdf2 %>%
    mutate(
      !!paste0(pollutant, "b") := cut(
        .data[[pollutant]], 
        breaks = quantile(.data[[pollutant]], probs = 0:2 / 2, na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2")
      )
    )
}


assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants2) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_htn_years, status_htn) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_htn_years, status_htn) ~", pollutant, "+ age + as.factor(sex)  + (1|id) +
                                                             logunits + logbmi + logsmok + 
                                                             simd2009v2rank"),
      paste("Surv(tte_htn_years, status_htn) ~", pollutant, "+ age + as.factor(sex)  + (1|id)
                                                                + logunits + logbmi + logsmok + 
                                                                simd2009v2rank + as.factor(diabetes) +
                                                                 Total_cholesterol + HDL_cholesterol")
    )
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

   #initialise lists for models
   model_results <- list()
   assump_results <- list()
   aic_list <- list()
for (formula in formula_list) {
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>%
                                                   mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>%
                                     as.data.frame() %>% 
                                     mutate(model = model_name, pollutant = pollutant)
aic_list[[model_name]] <- extract_coxme_aic(model, pollutant, model_name)
   # Increment model index
      model_index <- model_index + 1     
}
assump_test[[pollutant]] <- bind_rows(assump_results) #
cox_results[[pollutant]] <- bind_rows(model_results) 
aic_results[[pollutant]] <- bind_rows(aic_list)
  }

head(cox_results)
head(assump_test)
head(aic_results)

cox_results_df <- bind_rows(cox_results)
assump_test_df <- bind_rows(assump_test)
assump_test_df$term <- rownames(assump_test_df)
aic_results_df <- bind_rows(aic_results)

# #save out
fwrite(cox_results_df, file = "/Pollution_PCs/htn_coxme_bmipollcat.csv")
fwrite(assump_test_df, file = "/Pollution_PCs/htn_cox_assump_bmipollcat.csv")
fwrite(aic_results_df, file = "/Pollution_PCs/htn_cox_aic_bmipollcat.csv")



#_______________new coxdf for dementia analysis_____________#

#Load cox df 
cox_df <- fread("/Dementia/cox_df_sc_imp_16052025.csv")
table(cox_df$status)
cox_df$id <- as.character(cox_df$id)

#add PCs
coxdf2 <- left_join(cox_df, pc, by = "id")

#create categorical pollutant variable for each pollutant
coxdf2 <- coxdf
for (pollutant in pollutants1) {
  coxdf2 <- coxdf2 %>%
    mutate(
      !!paste0(pollutant, "b") := cut(
        .data[[pollutant]], 
        breaks = quantile(.data[[pollutant]], probs = 0:2 / 2, na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2")
      )
    )
}

assump_test <- list()
cox_results <- list()
aic_results <- list()

for (pollutant in pollutants2) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
     paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                        logunits + logbmi + logsmok + 
                                                        simd2009v2rank"), 
    paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id) +
                                                      logunits + logbmi + logsmok + simd2009v2rank + 
                                                      as.factor(e4_count)") 
    )
  # Convert formulas to formula objects
    formula_list <- lapply(formula_list, as.formula)
    
    # Model indexing
    model_index <- 1

   #initialise lists for models
   model_results <- list()
   assump_results <- list()
   aic_list <- list()
   model_summaries <- list()
   
for (formula in formula_list) {
model <- coxme(formula, varlist = kin_model*2, data = coxdf2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
# model_summaries[[pollutant]] <- summary(model)
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>% 
                                                mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>%
                                           as.data.frame() %>% 
                                           mutate(model = model_name, pollutant = pollutant)
aic_list[[model_name]] <- extract_coxme_aic(model, pollutant, model_name)
   # Increment model index
      model_index <- model_index + 1     
}
assump_test[[pollutant]] <- bind_rows(assump_results)
cox_results[[pollutant]] <- bind_rows(model_results) 
aic_results[[pollutant]] <- bind_rows(aic_list)
  }

head(cox_results)
head(assump_test)
head(aic_results)

cox_results_df <- bind_rows(cox_results)
assump_test_df <- bind_rows(assump_test)
assump_test_df <- assump_test_df %>% mutate(term = rownames(assump_test_df))
assump_test_df$term <- gsub("\\...\\d+", "", assump_test_df$term)
aic_results_df <- bind_rows(aic_results)

# #save out
fwrite(cox_results_df, file = "/Pollution_PCs/dementia_coxme_bmipollcat.csv")
fwrite(assump_test_df, file = "/Pollution_PCs/dementia_assump_bmipollcat.csv")
fwrite(aic_results_df, file = "/Pollution_PCs/dementia_cox_aic_bmipollcat.csv")

#models all run
##Load results
cvd_df <- fread("/Pollution_PCs/coxme_agebmipollcat_allcardiac.csv")
cvd_df <- cvd_df %>% mutate(outcome = "cvd")
htn_df <- fread("Pollution_PCs/htn_coxme_bmipollcat.csv")
htn_df <- htn_df %>% mutate(outcome = "htn")
mi_df  <- fread("/Pollution_PCs/coxme_agebmipollcat_mi_0707.csv")
mi_df <- mi_df %>% mutate(outcome = "mi")
ihd_df <- fread("/Pollution_PCs/coxme_agebmipollcat_ihd_0707.csv")
ihd_df <- ihd_df %>% mutate(outcome = "ihd")
stroke_df <- fread("/Pollution_PCs/coxme_agebmipollcat_stroke_0707.csv")
stroke_df <- stroke_df %>% mutate(outcome = "stroke")
dementia_df <- fread("/Pollution_PCs/dementia_coxme_bmipollcat.csv")
dementia_df <- dementia_df %>% mutate(outcome = "dementia")



#Filter to pollutant terms only
df <- rbind(cvd_df, htn_df, mi_df, ihd_df, stroke_df, dementia_df)
df <- df %>% 
  mutate(significance = ifelse(p.value < 0.05 & p.value > 0.05/2, "nominal", 
                          ifelse(p.value < 0.05/2, "bonf", "none")))

df <- df %>% filter(str_detect(term, pollutant))
savdf <- df %>%
  filter(model == "model_3") %>%
  mutate(
  pollutant = str_sub(pollutant, end = -2),
  term = pollutant)
savdf <- savdf %>% select(outcome, model, pollutant, everything())  
fwrite(savdf, file = "/Pollution_PCs/pc_results_forsupp.csv")

sig <- df %>% 
  group_by(outcome) %>% 
  dplyr::filter(significance == "bonf" & model == "model_3") %>%
  summarise(
    outcome = outcome,
    pollutant = str_sub(pollutant, end = -2),
    HR = estimate,
    CI = paste0(round(conf.low, 2), "-", round(conf.high, 2)),
    P_value = round(p.value, 3))


#filter to pollutant term only
df <- df %>%
  filter(term == paste0(pollutant, "Q2")) %>%
  mutate(significance = ifelse(p.value < 0.05/2, "P<0.025",
                          ifelse(p.value < 0.05 & p.value > 0.05/2, "P<0.05", "P>=0.05")))

df$pollutant <- str_sub(df$pollutant, end = -2)

library(paletteer)
#Plot with all models
p1 <- df3 %>%
ggplot(aes(y=pollutant, x=estimate, xmin=conf.low, xmax=conf.high, colour = model, shape = significance)) +
  geom_point(position = position_dodge(width = 0.75), size = 2) + 
  geom_errorbar(position = position_dodge(width = 0.75), height=.1) +
  scale_color_paletteer_d("PNWColors::Sunset2") +
  scale_shape_manual(values = c(17, 18, 16, 1)) +
    coord_cartesian(xlim = c(0.6, 2.6)) +
  scale_x_continuous(breaks = seq(0.6, 2.6, by = 0.4)) +
  labs(title='',
     x='Hazard Ratio [95% Confidence Interval]', 
     y = 'Outcome', colour = "model") +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_bw() + 
  theme(legend.position = "right",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title = element_text(size = 10),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(vjust = -0.5),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        strip.placement = "outside") +
  facet_wrap(~outcome,
  labeller = labeller(
    outcome = c("htn" = "Hypertension",
                "cvd" = "CVD",
                "stroke" = "Ischaemic Stroke",
                "mi" = "Myocardial Infarction",
                "ihd" = "Coronary Heart Disease",
                "dementia" = "Dementia")
  ))
ggsave(p1, file ="/Pollution_PCs/allmodels_summaryplot.pdf", width = 25, height = 30, units = "cm")

#Plot with final model only
# df$significance <- factor(df$significance, levels = c("none", "nominal", "yes"))
df$outcome <- factor(df$outcome, levels = c("dementia", "stroke", "htn", "ihd", "mi", "cvd"))

df$pollutant  <- fct_reorder(df$pollutant, df$estimate)

#Pick some colours
library(paletteer)

scale_colour_paletteer_d("MetBrewer::Benedictus")
scale_color_paletteer_d("MetBrewer::Benedictus")
scale_fill_paletteer_d("MetBrewer::Benedictus")
colour_list <- paletteer_d("MetBrewer::Benedictus")

col <- paletteer_d("NatParksPalettes::Glacier")
col <- col[1:4]

#add column with case numbers
df <- df %>% 
  mutate(case_lable = ifelse(outcome == "cvd", "Cases:905 \nControls:12118",
                        ifelse(outcome == "htn", "Cases:1689 \nControls:19064",
                          ifelse(outcome == "mi", "Cases:145 \nControls:12878",
                            ifelse(outcome == "ihd", "Cases:530 \nControls:12493",
                              ifelse(outcome == "stroke", "Cases:159 \nControls:12864", "Cases:323 \nControls:10416"))))))

# Decide where to place the label inside each facet:
# - x position just inside the left x-limit
# - y position at the top (or bottom) pollutant, nudged a bit to avoid overlap
x_left <- 2.2  # inside your coord_cartesian(xlim = c(0.6, 2.6))
poll_levels <- levels(fct_inorder(df$pollutant))
# y_top <- poll_levels[length(poll_levels)]      # top row
y_bottom <- poll_levels[1]                   # or bottom row if you prefer

# One row per facet (outcome) with a single case_lable
labels_df <- df %>%
  group_by(outcome) %>%
  summarise(case_lable = first(na.omit(case_lable)), .groups = "drop") %>%
  mutate(
    estimate  = x_left,
    pollutant = y_bottom,
  )



p2 <- df %>% filter(model == "model_3") %>%
  ggplot(aes(y=pollutant, x=estimate, xmin=conf.low, xmax=conf.high, shape = significance, colour = significance)) +
  geom_point(position = position_dodge(width = 0.75), size = 2.75) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.2, linewidth = 1.05) +
  scale_color_manual(name = "Significance", values = col) +
  scale_shape_manual(name = "Significance", values = c(17, 18, 16, 1)) +
    coord_cartesian(xlim = c(0.4, 2.2)) +
  scale_x_continuous(breaks = seq(0.4, 2.2, by = 0.4)) +
    labs(title='',
     x='Hazard Ratio [95% Confidence Interval]', 
     y = 'Outcome', colour = "Significance") +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_bw() + 
  theme(legend.position = "right",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title = element_text(size = 14),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        strip.placement = "outside") +
  facet_wrap(~outcome,
  labeller = labeller(
    outcome = c("htn" = "Hypertension",
                "cvd" = "CVD",
                "stroke" = "Ischaemic Stroke",
                "mi" = "Myocardial Infarction",
                "ihd" = "Coronary Heart Disease",
                "dementia" = "Dementia")
  )) +
    geom_text(
    data = labels_df,
    aes(x = estimate, y = pollutant, label = case_lable),
    inherit.aes = FALSE,
    hjust = 0.5,              # left-align at x_left
    vjust = 9,           # nudge a little above the chosen row
    size = 3
  ) +
  theme(strip.text.x = element_text(size = 16))

ggsave(p2, file = "/Pollution_PCs/model3_summaryplot2.pdf", width = 40, height = 30, units = "cm")


#Also save out df of proportional hazards for supplementary
strokeassump <- fread("/Pollution_PCs/cox_assump_agebmipollcat__stroke_0707.csv")
strokeassump <- strokeassump %>% mutate(outcome = "stroke")
cardiacassump <- fread("/Pollution_PCs/cox_assump_agebmipollcat_allcardiac.csv")
cardiacassump <- cardiacassump %>% mutate(outcome = "cardiac")
ihdassump <- fread("/Pollution_PCs/cox_assump_agebmipollcat_ihd_0707.csv")
ihdassump <- ihdassump %>% mutate(outcome = "ihd")
miassump <- fread("/Pollution_PCs/cox_assump_agebmipollcat_mi_0707.csv")
miassump <- miassump %>% mutate(outcome = "mi")
demassump <- fread("/Pollution_PCs/dementia_assump_bmipollcat.csv")
demassump <- demassump %>% mutate(outcome = "dementia")
htnassump <- fread("/Pollution_PCs/htn_cox_assump_bmipollcat.csv")
htnassump <- htnassump %>% mutate(outcome = "htn")

assumpdf <- rbind(cardiacassump, ihdassump, miassump, strokeassump, htnassump, demassump)
assumpdf <- assumpdf %>% filter(model == "model_3") %>%
  rename(PC = pollutant) %>%
  mutate(PC = str_sub(PC, end = -2),
         term = gsub("\\...\\d+", "", term))
assumpdf <- assumpdf %>% filter(term == "GLOBAL" | str_detect(term, PC))
assumpdf <- assumpdf %>% select(outcome, PC, term, model, everything())
assumpdf %>% filter(p < 0.05)
fwrite(assumpdf, file = "/Pollution_PCs/PCs_health_coxassump.csv")
