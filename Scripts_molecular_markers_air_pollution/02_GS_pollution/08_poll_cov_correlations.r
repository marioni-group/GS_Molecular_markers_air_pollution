#pollution x covariate correlations

# Exploratory SIMD analysis
library(tidyverse)
library(data.table)
library(here)
library(Hmisc)
library(ComplexHeatmap)
library(here)
library(broom)

#4. Plot pollutant x SIMD

#Load SIMD

simd <- fread("2024-10-18_simd09_all_categories.csv") 

#Load Pollutant - annual data
annual_df <- fread("annual_avg_exp_foranalysis_scaled.csv")
head(annual_df)

#load simd data and merge
simd <- simd %>% rename(ID = id)
annual_df3 <- annual_df %>%
  dplyr::filter(ID %in% simd$ID) 

annual_exp_simd <- left_join(annual_df3, simd, by = "ID")
annual_exp_simd %>% group_by(ID) %>% summary()
annual_exp_simd <- annual_exp_simd %>% drop_na() 

#Load covariate data and merge
#Try with covariates
cov <- fread("covariates_knnimpute.csv") 

#add in year of appointment
extra <- fread("GS_appt_dates.txt")
extra <- extra %>% dplyr::select(id, ym) %>%
  dplyr::mutate(y_appt = ym(ym)) %>%
  dplyr::mutate(y_appt = year(y_appt))
head(extra)

#set-up labels for covariates
colnames(simd)
simd_comp <- c("age", "logsmok", "logbmi", "logunits", "y_appt", "simd2009v2rank", "simd2009v2_inc_rank",  "simd2009v2_emp_rank", "simd2009v2_hlth_rank", "simd2009v2_educ_rank","simd2009v2_house_rank", "simd2009v2_access_rank", "simd2009v2_crime_rank")
simd_labels <- c("Age", "Smoking", "BMI", "Alcohol", "SIMD rank", "year_range", "Income rank", "Employment rank", "Health rank", "Education rank", "House rank", "Access rank", "Crime rank")
head(annual_df3)

#select simd aspects of interest
simd_select <- simd %>% dplyr::select(any_of(c("ID", simd_comp)))  

#select covariates of interest & set-up for joining
cov_select <- cov %>% dplyr::select(any_of(c("id", simd_comp))) 
cov_select <- cov_select %>% dplyr::rename(ID = id)
cov_select <- cov_select %>% dplyr::select(-simd2009v2rank)
simd_select$ID <- as.character(simd_select$ID)
annual_df3$ID = as.character(annual_df3$ID)
cov_select$ID <- as.character(cov_select$ID)

#Join
simd_poll_corr <- left_join(annual_df3, cov_select, by = "ID")
simd_poll_corr <- left_join(simd_poll_corr, simd_select, by = "ID")
extra$ID <- as.character(extra$id)
simd_poll_corr <- left_join(simd_poll_corr, extra, by = "ID")
simd_poll_corr2 <- simd_poll_corr %>% dplyr::select(-c(ID, id, ym))

#Correlate Pollutant & SIMD variables & other covariates
#Run correlation
corroutput <- Hmisc::rcorr(as.matrix(simd_poll_corr2), type = "pearson") #year range not included
rmat <- corroutput$r
pmat <- corroutput$P

colnames(rmat)
rownames(rmat)
heatmap_labels <- c("PM25", "NO2", "NO", "O3", "PM10", 
                    "SO2", "NO3_C","NO3_F", "age", "smoking",
                     "BMI", "alcohol", "SIMD rank", "Income rank", 
                     "Employment rank", "Health rank", "Education rank", 
                     "House rank", "Access rank", "Crime rank", "Appt year")


rmat2 <- as.matrix(rmat)
rmat2[upper.tri(rmat)] <- NA

png("pollutant_cov_simd_htmp_18556impt.png")
Heatmap(rmat2,
    cell_fun = function(j, i, x, y, width, height, fill) {
        if(!is.na(rmat2[i,j])) {
        grid.text(sprintf("%.1f", rmat2[i, j]), x, y, gp = gpar(fontsize = 10))}},
  na_col = "white",
  name = "Correlation",
  row_order = rownames(rmat),
  column_order = colnames(rmat),
  row_labels = heatmap_labels,
  column_labels = heatmap_labels,
  column_title = "",
  row_title = "",
  row_names_side = "left",
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    title = "Correlation",
    direction = "vertical",
    just = c("top", "left"))
    )
dev.off()
