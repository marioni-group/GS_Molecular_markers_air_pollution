#Run cox - models for all-cause dementia
library(tidyverse)
library(data.table)
library(coxme)

#Load cox df 
cox_df <- fread("Dementia/cox_df_sc_imp_08072025.csv")

#load kinship matrix
kin_model <- readRDS("Dementia/kin_model_0901.rds")


#function to extract coxme results 
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


#Set up variables to iterate over
pollutants <- colnames(cox_df)[18:25]
pollutants


#Create high vs low pollutants 
cox_df2 <- cox_df
for (pollutant in pollutants) {
  cox_df2 <- cox_df2 %>%
    mutate(
      !!paste0(pollutant, "d") := cut(
        .data[[pollutant]], 
        breaks = quantile(.data[[pollutant]], probs = 0:2 / 2, na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2")
      )
    )
}

pollutants1 <- c("PM25d", "NO2d",  "NOd",  "O3d",  "PM10d", "SO2d", "NO3_Cd", "NO3_Fd")




assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors
#edit models for follow up years as required
for (pollutant in pollutants1) {
    print(pollutant)
    #model list
   formula_list <- list(
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + logunits + logbmi + logsmok + simd2009v2rank + as.factor(e4_count)"), 
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id) + logunits + logbmi + logsmok + simd2009v2rank +  as.factor(e4_count) + Total_cholesterol + HDL_cholesterol + as.factor(diabetes) + as.factor(hypertension)")
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
model <- coxme(formula, varlist = kin_model*2, data = cox_df)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% as.data.frame() %>% mutate(model = model_name, pollutant = pollutant)
aic_list[[model_name]] <- extract_coxme_aic(model, pollutant, model_name)
   # Increment model index
      model_index <- model_index + 1     
}
assump_test[[pollutant]] <- bind_rows(assump_results)
cox_results[[pollutant]] <- bind_rows(model_results, .id = "model") 
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

#save out
fwrite(cox_results_df, file = "Dementia/cox_results_pollcat_full_0807.csv")
fwrite(assump_test_df, file = "Dementia/cox_assump_pollcat_full_0807.csv")
fwrite(aic_results_df, file = "Dementia/cox_aic_pollcat_full_0807.csv")

