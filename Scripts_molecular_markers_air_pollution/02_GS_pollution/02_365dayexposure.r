#Pollution summaries using 2005 -2011 data to generate 365 day data

library(tidyverse)
library(data.table) 
library(lubridate) 

#Load data
#Load Data: saved in two chunks currently
results1 <- readRDS("/EMEP4uk_results_part1.rds")
results2 <- readRDS("/EMEP4uk_results_part2.rds")

#create dataframe of results
df1 <- dplyr::bind_rows(results1)
df2 <- dplyr::bind_rows(results2)

#combine together
full_results <- rbind(df1, df2)
gc()

#extract 365 day data for relevant time periods for each ID
#use GS appointment dates as end date for pollution capture e.g. 
app_dates <- fread("GS_appt.txt")
app_dates$appt <- gsub(" .*", "", app_dates$appt) #remove time from appt as don't require this
head(app_dates)


# Create an interval of 365 days
time_pd <- days(365)
app_dates$appt <- as.Date(app_dates$appt)
app_dates1 <- app_dates %>%
  mutate(start_date = appt - time_pd,
         id = as.character(id))
head(app_dates1)
sum(is.na(app_dates1))
app_dates1 <- app_dates1 %>%
  rename(ID = id) %>%
  mutate(ID = as.integer(ID))
  dim(app_dates1) #24081

#add in dates to results df
full_results <- left_join(full_results, app_dates1, by = "ID")

#check methylation IDs
meth_ids <- fread("2024-10-31_meth_ids.csv")
#Load list of IDs to keep
gs_ids <- read.table("/GS_ids_20241007.txt") 
full_resultsids <- unique(full_results$ID)

#Filter to remove ids not in gs_ids
full_results <- full_results %>%
  dplyr::filter(ID %in% gs_ids$V1)
length(unique(full_results$ID))

#subset data to relevant dates for each ID
full_results <- full_results %>%
  rename(measurement_day = Date)

full_results <- full_results %>%
  mutate(measurement_day = as.Date(measurement_day))

#_________________________________________________________________________
#filter to 365 days before the appointment (not including appointment day)
library(foreach)
library(doParallel)

# Convert to data.table
setDT(full_results)

# Set up cluster
cl <- makeCluster(12)
registerDoParallel(cl)

# Parallel processing with foreach
filtered_results_list <- foreach(split_data = iter(split(full_results, by = "ID")), .packages = "data.table") %dopar% {
  split_data[measurement_day >= start_date & measurement_day < appt]
}

# Stop cluster
stopCluster(cl)

# Combine the results
filtered_results <- rbindlist(filtered_results_list)

#visually inspect a few IDs
filtered_results[ID == "33", .N, by = pollutant] #365 entries for each pollutant 

# Count NAs for each column
na_counts <- filtered_results[, lapply(.SD, function(x) sum(is.na(x)))]
 #108040 measurements missing
na_counts_full <-  full_results[, lapply(.SD, function(x) sum(is.na(x)))]

#Save out 365day full results
saveRDS(filtered_results, file = "/Pollution/EMEP4uk_gs/daily_365_pollution.RDS")  

#_________________________________________________________________________________________
#Create summary data

#add year column
filtered_results <- filtered_results %>%
 mutate(Year = year(measurement_day))

#Create 365day average data
annual_df <- filtered_results[, .(
  days_365_avg = mean(Measurement, na.rm = TRUE),
  standard_dev = sd(Measurement, na.rm = TRUE),
  year_range = list(unique(Year))
), by = .(ID, pollutant)]

#________________________________________________________________________________________________
#Clean up pollution exposure df
#Complete pollutant abbreviations
annual_df <- annual_df %>% dplyr::rename(raster_entry = pollutant)
annual_df$Pollutant <- annual_df$raster_entry
head(annual_df)
annual_df2 <- annual_df %>%
  mutate(Pollutant = ifelse(raster_entry == "SURF_ug_NO", "NO", 
                      ifelse(raster_entry == "SURF_ug_NO3_C", "NO3_C", 
                        ifelse(raster_entry == "SURF_ug_NO3_F", "NO3_F",
                          ifelse(raster_entry == "SURF_ug_PM25_rh50", "PM25",
                            ifelse(raster_entry == "SURF_ppb_O3", "O3", 
                      ifelse(raster_entry == "SURF_ug_NO2", "NO2", 
                        ifelse(raster_entry == "SURF_ug_PM10_rh50", "PM10",
                            ifelse(raster_entry == "SURF_ug_SO2", "SO2", Pollutant)))))))))

#check all years present
unique(annual_df2$year_range) 

#check NA values
sum(is.na(annual_df2))
lapply(annual_df2, function(x) sum(is.na(x))) 
annual_df3 <- annual_df2 %>% drop_na() #22071 

#Save out scaled
annual_df4 <- annual_df3 %>%
  dplyr::select(ID, Pollutant, days_365_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )
scaled_pollution <- annual_df4  
scaled_pollution[2:9] <- scale(annual_df4[2:9])   


fwrite(annual_df3, file = "/Pollution/EMEP4uk_gs/annual_avg_exp_foranalysis_ntr.csv")
fwrite(scaled_pollution, file = "/Pollution/EMEP4uk_gs/annual_avg_exp_foranalysis_scaled.csv")


####Correlations
#Calculate between-pollutant correlations
df <- fread("/Pollution/EMEP4uk_gs/annual_avg_exp_foranalysis_ntr.csv")
df2 <- df %>% dplyr::select(ID, days_365_avg, Pollutant)
head(df2)

df2 %>%
  dplyr::summarise(n = dplyr::n(), .by = c(Pollutant)) %>%
  dplyr::filter(n > 1L) 

df3 <- df2 %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )

cor_poll <- cor(df3[2:9])
pdf("/Pollution/EMEP4uk_gs/between_poll_htmp_365avg.pdf")
poll_heatmap <- ComplexHeatmap::Heatmap(cor_poll)
ComplexHeatmap::draw(poll_heatmap)
dev.off()

fwrite(cor_poll, file = "/Pollution/EMEP4uk_gs/between_poll_corr_matrix.csv")


#______________________________________________________________________******for separate script********
#add in 6 month and 7 year data  
sevenyr <- fread("Pollution/EMEP4uk_gs/sevenyr_foranalysis_notsc.csv")
sevenyr <- sevenyr %>% select(ID, avg_7yrs, Pollutant)
sixmnth <- fread("Pollution/EMEP4uk_gs/6month_avg_exp_foranalysis_ntr.csv")
sixmnth <- sixmnth %>% pivot_longer(!ID, names_to = "Pollutant", values_to = "avg_6mnth")
annual_df3 <- annual_df2 %>% select(ID, Pollutant, days_365_avg)

df <- left_join(sixmnth, annual_df3, by = c("ID", "Pollutant"))
df <- left_join(df, sevenyr, by = c("ID", "Pollutant"))
df2 <- df %>% pivot_longer(!c(ID, Pollutant), names_to = "duration", values_to = "exposure")

df2$duration <- factor(df2$duration, levels = c("avg_6mnth", "days_365_avg", "avg_7yrs"))
str(df2)

p <- df2 %>%
  ggplot(aes(x = duration, y = exposure, fill = Pollutant)) +
  geom_violin(trim = TRUE, draw_quantiles =  c(0.25, 0.5, 0.75), alpha = 0.8) +
  scale_fill_paletteer_d("MetBrewer::Tam") +
  scale_x_discrete(labels = c("6 month", "1 year", "7 years")) +
  # geom_jitter(colour = "black", size = 0.2, alpha = 0.3) +
  theme_minimal() +
  labs(y = "Average exposure over varying time frames", x = "Duration") +
  facet_wrap(~Pollutant, scales = "free_y", nrow = 2) +
  theme(legend.position = "bottom")
ggsave(p, file = "Pollution/EMEP4uk_gs/exp_summaries_vplots.pdf", width = 30, height = 15, units = "cm")


#Create a summary table
library(arsenal)
  
annual_df2 <- annual_df2 %>%
  mutate(y_appt = gsub("*....\\|", "", year_range))  
unique(annual_df2$year_range)


control <- tableby.control(,    # Use ANOVA for numeric tests
  total = TRUE,              # Include a total column
  numerics = list(
    mean = FALSE,
    Range = FALSE,
    median = TRUE,           # Show median
    q1q3 = TRUE)) 

tab1 <- summary(tableby(y_appt ~ days_365_avg, 
                          data = annual_df2, 
                          strata = Pollutant, 
                          control = control), text = TRUE)
fwrite(tab1, file = "Pollution/EMEP4uk_gs/summary_table.csv")


summary_data <- annual_df2 %>%
  group_by(y_appt, Pollutant) %>%
  summarise(
    count = n(),
    median_value = round(median(days_365_avg, na.rm = TRUE), 2),
    IQR_Value = round(IQR(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(median_value, "[", IQR_Value, "]")
  )
summary_by_pollutant <- annual_df2 %>%
  group_by(Pollutant) %>%
  summarise(
    count = n(),
    median_value = round(median(days_365_avg, na.rm = TRUE), 2),
    IQR_Value = round(IQR(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(median_value, "[", IQR_Value, "]")
  )

summary_wide <- summary_data %>%
  select(-median_value, -IQR_Value, -count) %>%
  pivot_wider(names_from = y_appt, values_from = summary_text, values_fill = "")

fwrite(summary_wide, file = "Pollution/EMEP4uk_gs/summary_table2.csv" )
sum_dat <- fread("Pollution/EMEP4uk_gs/summary_table2.csv")

#Also create mean/S.d. summary table
library(arsenal)
  
annual_df2 <- annual_df2 %>%
  mutate(y_appt = gsub("*....\\|", "", year_range))  
unique(annual_df2$year_range)


control <- tableby.control(,    # Use ANOVA for numeric tests
  total = TRUE,              # Include a total column
  numerics = list(
    meansd = TRUE)) 

tab1 <- summary(tableby(y_appt ~ days_365_avg, data = annual_df2, strata = Pollutant, control = control), text = TRUE)
fwrite(tab1, file = "Pollution/EMEP4uk_gs/summary_table_meansd.csv")


summary_data <- annual_df2 %>%
  group_by(y_appt, Pollutant) %>%
  summarise(
    count = n(),
    mean_value = round(mean(days_365_avg, na.rm = TRUE), 2),
   sd_value = round(sd(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(mean_value, "[", sd_value, "]")
  )
summary_by_pollutant <- annual_df2 %>%
  group_by(Pollutant) %>%
  summarise(
    count = n(),
    mean_value = round(mean(days_365_avg, na.rm = TRUE), 2),
    sd_Value = round(sd(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(mean_value, "[", sd_Value, "]")
  )

fwrite(summary_by_pollutant, file = "Pollution/EMEP4uk_gs/pollutant_summary_meansd.csv")

