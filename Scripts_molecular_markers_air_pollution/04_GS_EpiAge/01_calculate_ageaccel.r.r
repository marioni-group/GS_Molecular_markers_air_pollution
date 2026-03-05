#Calculating epi-age acceleration

# model <- Epi-age ~ age
# Age-acceleration = residuals(model)
# To look at: GrimAge v2, PhenoAge, DunedinPACE

clock_data <- read.csv("gs_clocks.csv")
target <- readRDS("GS20k_Targets_18869.rds")
pheno <- readRDS("GS_phenos_internal_23Aug2023_REM.rds")
target_age <- target %>% dplyr::select(Sample_Sentrix_ID, age)

clock2 <- clock_data %>% dplyr::select(X.1, PhenoAge, DunedinPACE, DNAmGrimAge.1) %>% #GrimAge.1 = GrimAgev2
  rename(Sample_Sentrix_ID = X.1)
clock2 <- left_join(target_age, clock2, by = "Sample_Sentrix_ID")  
dim(clock2) 

clocks <- c("PhenoAge", "DunedinPACE", "DNAmGrimAge.1")

clock_resid <- clock2
clock_resid <- clock_resid %>% drop_na()
clock2 <- clock2 %>% drop_na()
ageaccel_models <- list()

for (clock in clocks) {
 print(clock)
 model <- lm(clock2[[clock]] ~ age, data = clock2)
 clock_resid[[clock]] <- paste(resid(lm(clock2[[clock]] ~ age, data = clock2)))
 ageaccel_models[[clock]] <- broom::tidy(model)
}

#convert character output from loop to numeric
clock_resid <- clock_resid %>%
  purrr::modify_at(c(3:9), as.numeric)
clock_resid2 <- clock_resid %>%
  dplyr::rename_with(~paste0(.x, '_accel'), 3:9)

fwrite(clock_resid2, file = "age_acceleration_clocks.csv")