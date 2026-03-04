#03_htn_coxmodels_run_070725

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


#Updated: load files from pre-prepared scripts
library(tidyverse)
library(data.table)
library(coxme)
# #load kinship matrix
kin_model <- readRDS("kin_model_0901.rds")

#Load cox df 

cox_df <- fread("htn_cox_df_bmifilt_sc_imp_070725.csv")

dim(cox_df)
table(cox_df$status_htn) #1689
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


assump_test <- list()
cox_results <- list()
aic_results <- list()

# Iterate over pollutants and predictors

for (pollutant in pollutants1) {
    print(pollutant)
    #model list
    formula_list <- list(
      paste("Surv(tte_htn_years, status_htn) ~", pollutant, "+ age + as.factor(sex) + (1|id)"),
      paste("Surv(tte_htn_years, status_htn) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                              logunits + logbmi + logsmok + simd2009v2rank"),
      paste("Surv(tte_htn_years, status_htn) ~", pollutant, "+ age + as.factor(sex)  + (1|id) + 
                                                              logunits + logbmi + logsmok + simd2009v2rank +
                                                               as.factor(diabetes) + Total_cholesterol + HDL_cholesterol")
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

head(cox_results)
head(assump_test)
head(aic_results)

cox_results_df <- bind_rows(cox_results)
assump_test_df <- bind_rows(assump_test)
assump_test_df$term <- rownames(assump_test_df)
aic_results_df <- bind_rows(aic_results)

# #save out
fwrite(cox_results_df, file = "htn_coxme_bmipollcat_0807.csv")
fwrite(assump_test_df, file = "htn_cox_assump_bmipollcat_0807.csv")
fwrite(aic_results_df, file = "htn_cox_aic_bmipollcat_0807.csv")



#Check proportional hazards in coxph
##Explore time-varying effects
# Collect results
hr_all <- map_dfr(pollutants1, function(pollutant) {
  map_dfr(5:18, function(y) {
  
  # Truncate time and status
  df_temp <- cox_df2 %>%
    mutate(
      tte_temp = pmin(tte_htn_years, y),
      status_temp = ifelse(tte_htn_years > y, 0, status_htn)
    )
  
  # Fit Cox model
formula_text <- paste0("Surv(tte_temp, status_temp) ~", pollutant, "+ age + as.factor(sex)
                                                             + logunits + logbmi + logsmok + 
                                                             simd2009v2rank + as.factor(diabetes) +
                                                              Total_cholesterol + HDL_cholesterol")
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

fwrite(hr_all, file = "htn_coxph_bmipollcat_0807.csv")
head(hr_all)
table(hr_all$ph_status)
failed <- hr_all %>% filter(ph_status == "PH Failed")
table(failed$pollutant)
failed


#Plot
coxhr <- ggplot(hr_all, aes(x = year, y = estimate, color = ph_status)) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.4) +
  geom_text(aes(y = 0.2, label = n_label), angle = 90, size = 3, vjust = 0.5, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("PH Passed" = "blue", "PH Failed" = "red", "Unknown" = "gray")) +
  labs(
    title = "Hazard Ratios of HTN Over Time for Multiple Pollutants",
    x = "Years of Follow-up",
    y = "Hazard Ratio (HR)",
    color = "PH Assumption"
  ) +
  facet_wrap(~ pollutant, scales = "free_y") +
  theme_minimal()

ggsave(coxhr, file = "htn_coxph_bmipollcontHRplot_0807.pdf", height = 20, width = 30, units = "cm")
