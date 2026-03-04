#Updated Pollution_cvd cox models

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
library(coxme)

# #load kinship matrix
kin_model <- readRDS("kin_model_0901.rds")

#Load cox df 

filename <- "cvd_cox_df_agebmifilt_sc_imp_070725.csv" #edit for relevant input
path <- ""
cox_df <- fread(paste0(path, filename))


dim(cox_df)
table(cox_df$status)
summary(cox_df$tte_years)

pollutants <- c("PM25", "NO2",  "NO",  "O3",  "PM10", "SO2", "NO3_C", "NO3_F")

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


##Explore time-varying effects
# Collect results
hr_all <- map_dfr(pollutants, function(pollutant) {
  map_dfr(5:18, function(y) {
  
  # Truncate time and status
  df_temp <- cox_df2 %>%
    mutate(
      tte_temp = pmin(tte_years, y),
      status_temp = ifelse(tte_years > y, 0, status)
    )
  
  # Fit Cox model
formula_text <- paste0("Surv(tte_temp, status_temp) ~", pollutant, " + age + as.factor(sex) + logunits + logbmi + logsmok + simd2009v2rank + as.factor(diabetes) + as.factor(hypertension) + Total_cholesterol + HDL_cholesterol")
formula_obj <- as.formula(formula_text)

#fit model
model <- coxph(formula_obj, data = df_temp)

    # PH assumption test
ph_p <- tryCatch({
      ph <- cox.zph(model)
      row_match <- grep(paste0("^", pollutant), rownames(ph$table))
      if (length(row_match) == 0) NA else max(ph$table[row_match, "p"], na.rm = TRUE)
    }, error = function(e) NA)
    
    ph_status <- ifelse(is.na(ph_p), "Unknown", ifelse(ph_p >= 0.05, "PH Passed", "PH Failed"))

    # Count cases and controls at this time point
  n_cases <- sum(df_temp$status_temp == 1, na.rm = TRUE)
  n_controls <- sum(df_temp$status_temp == 0, na.rm = TRUE)
  n_text <- paste0("n=", n_cases, "/", n_controls)


  # Extract HR and CI for your_var
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl(paste0("^", pollutant), term)) %>%
    mutate(year = y,
    ph_status = ph_status,
    cases = n_cases,
    controls = n_controls,
    n_label = n_text,
    pollutant = pollutant)

})
})

fwrite(hr_all, file = "cvd_coxph_agebmipollcont_0707.csv")

#Plot
coxhr <- ggplot(hr_all, aes(x = year, y = estimate, color = ph_status)) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.4) +
  geom_text(aes(y = 0.2, label = n_label), angle = 90, size = 3, vjust = 0.5, hjust = 0.2, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("PH Passed" = "blue", "PH Failed" = "red", "Unknown" = "gray")) +
  labs(
    title = "Hazard Ratios of CVD Over Time for Multiple Pollutants",
    x = "Years of Follow-up",
    y = "Hazard Ratio (HR)",
    color = "PH Assumption"
  ) +
  facet_wrap(~ pollutant, scales = "free_y") +
  theme_minimal()

ggsave(coxhr, file = "coxph_agebmipollcontHRplot_020126.pdf", height = 30, width = 30, units = "cm")


#Repeat for categorical pollutant
##Explore time-varying effects
# Collect results
hr_all <- map_dfr(pollutants1, function(pollutant) {
  map_dfr(5:18, function(y) {
  
  # Truncate time and status
  df_temp <- cox_df2 %>%
    mutate(
      tte_temp = pmin(tte_years, y),
      status_temp = ifelse(tte_years > y, 0, status)
    )
  
  # Fit Cox model
formula_text <- paste0("Surv(tte_temp, status_temp) ~", pollutant, 
                        " + age + as.factor(sex) + logunits + logbmi +
                        logsmok + simd2009v2rank + as.factor(diabetes) + 
                        as.factor(hypertension) + Total_cholesterol + HDL_cholesterol")
formula_obj <- as.formula(formula_text)

#fit model
model <- coxph(formula_obj, data = df_temp)

    # PH assumption test
ph_p <- tryCatch({
      ph <- cox.zph(model)
      row_match <- grep(paste0("^", pollutant), rownames(ph$table))
      if (length(row_match) == 0) NA else max(ph$table[row_match, "p"], na.rm = TRUE)
    }, error = function(e) NA)
    
    ph_status <- ifelse(is.na(ph_p), "Unknown", ifelse(ph_p >= 0.05, "PH Passed", "PH Failed"))

    # Count cases and controls at this time point
  n_cases <- sum(df_temp$status_temp == 1, na.rm = TRUE)
  n_controls <- sum(df_temp$status_temp == 0, na.rm = TRUE)
  n_text <- paste0("n=", n_cases, "/", n_controls)


  # Extract HR and CI for your_var
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl(paste0("^", pollutant), term)) %>%
    mutate(year = y,
    ph_status = ph_status,
    cases = n_cases,
    controls = n_controls,
    n_label = n_text,
    pollutant = pollutant)

})
})


#Plot
coxhr <- ggplot(hr_all, aes(x = year, y = estimate, color = ph_status)) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.4) +
  geom_text(aes(y = 0.2, label = n_label), angle = 90, size = 3, vjust = 0.5, hjust = 0.2, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("PH Passed" = "blue", "PH Failed" = "red", "Unknown" = "gray")) +
  labs(
    title = "Hazard Ratios of CVD Over Time for Multiple Pollutants (high/low)",
    x = "Years of Follow-up",
    y = "Hazard Ratio (HR)",
    color = "PH Assumption"
  ) +
  facet_wrap(~ factor(pollutant, levels = c("NOd", "NO2d", "NO3_Cd", "NO3_Fd", "O3d", "PM10d", "PM25d", "SO2d")), scales = "free_y") +
  theme_minimal()

ggsave(coxhr, file = "coxph_agebmipollcatHRplot_020126.pdf", height = 30, width = 30, units = "cm")


# #Based on above exploration -> use pollutants divided into high vs. low, age & bmi filtered data
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


# #Set up variables to iterate over
pollutants1 #pollutants stratified into high vs. low


assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants1) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
              logunits + logbmi + logsmok + simd2009v2rank"),
      paste("Surv(tte_years, status) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
              logunits + logbmi + logsmok + simd2009v2rank + as.factor(diabetes) + 
              Total_cholesterol + HDL_cholesterol + as.factor(hypertension)")
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
model <- coxme(formula, varlist = kin_model*2, data = cox_df2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>%
                                mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% 
                                as.data.frame() %>% 
                                mutate(model = model_name, 
                                        pollutant = pollutant)
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
assump_test_df$term <- rownames(assump_test_df)
aic_results_df <- bind_rows(aic_results)

# #save out
fwrite(cox_results_df, file = "coxme_cvdagebmipollcat_0707.csv")
fwrite(assump_test_df, file = "coxme_assump_cvdagebmipollcat_0707.csv")
fwrite(aic_results_df, file = "coxme_aic_cvdagebmipollcat_0707.csv")


#__________________________
#Also run for sub-categories of CVD
#__________________________
#Ischaemic stroke
# #Set up variables to iterate over
pollutants1 <- c("PM25b", "NO2b", "NOb", "O3b", "PM10b", "SO2b", "NO3_Cb", "NO3_Fb")#pollutants stratified into high vs. low
pollutants <- c("PM25", "NO2", "NO", "O3", "PM10", "SO2", "NO3_C", "NO3_F")

#load cox_df
cox_df2 <- fread("cvd_cox_df_agebmifilt_sc_imp_070725.csv")

formula <- as.formula("Surv(tte_years, status_stroke) ~ PM25b + age + 
                        as.factor(sex)  + (1|id) + logunits + logbmi + 
                        logsmok + simd2009v2rank + as.factor(diabetes) + 
                        Total_cholesterol + HDL_cholesterol + 
                        as.factor(hypertension)")
model <- coxme(formula, varlist = kin_model*2, data = cox_df2)
summary(model)


#Add in categorised pollutants
for (pollutant in pollutants) {
  cox_df2 <- cox_df2 %>%
    mutate(
      !!paste0(pollutant, "b") := cut(
        .data[[pollutant]], 
        breaks = quantile(.data[[pollutant]], probs = 0:2 / 2, na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2")
      ))}

#make sure outcomes have NAs as 0 
cox_df2 <- cox_df2 %>%
  mutate(status_stroke = coalesce(status_stroke, 0),
         status_IHD = coalesce(status_IHD, 0),
         status_MI = coalesce(status_MI, 0))

assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants1) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status_stroke) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status_stroke) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                logunits + logbmi + logsmok + simd2009v2rank"),
      paste("Surv(tte_years, status_stroke) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                            logunits + logbmi + logsmok + simd2009v2rank +
                                                            as.factor(diabetes) + Total_cholesterol + 
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
model <- coxme(formula, varlist = kin_model*2, data = cox_df2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>%
                               mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% as.data.frame() %>% 
                                mutate(model = model_name, 
                                pollutant = pollutant)
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

fwrite(cox_results_df, file = "coxme_agebmipollcat_stroke_0707.csv")
fwrite(assump_test_df, file = "cox_assump_agebmipollcat__stroke_0707.csv")
fwrite(aic_results_df, file = "cox_aic_agebmipollcat_stroke_0707.csv")

#__________________________
#IHD
assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants1) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status_IHD) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status_IHD) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                          logunits + logbmi + logsmok + 
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
model <- coxme(formula, varlist = kin_model*2, data = cox_df2)
model_name <- paste("model", model_index, sep = "_")
model_output_name <- paste(pollutant, model_index, sep = "_")
model_results[[model_name]] <- extract_coxme_results(model, pollutant, model_name) %>%
                                 mutate(model = model_name)
assump_results[[model_name]] <- cox.zph(model)$table %>% as.data.frame() %>% 
                                                        mutate(model = model_name, 
                                                        pollutant = pollutant)
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

fwrite(cox_results_df, file = "coxme_agebmipollcat_ihd_0707.csv")
fwrite(assump_test_df, file = "cox_assump_agebmipollcat_ihd_0707.csv")
fwrite(aic_results_df, file = "cox_aic_agebmipollcat_ihd_0707.csv")

#__________________________
# #MI
assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants1) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_years, status_MI) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_years, status_MI) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                          logunits + logbmi + logsmok + simd2009v2rank"),
      paste("Surv(tte_years, status_MI) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                        logunits + logbmi + logsmok + simd2009v2rank + 
                                                        as.factor(diabetes) + Total_cholesterol + 
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
model <- coxme(formula, varlist = kin_model*2, data = cox_df2)
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

fwrite(cox_results_df, file = "coxme_agebmipollcat_mi_0707.csv")
fwrite(assump_test_df, file = "cox_assump_agebmipollcat_mi_0707.csv")
fwrite(aic_results_df, file = "ox_aic_agebmipollcat_mi_0707.csv")
