#Pollution summaries using 2005 -2011 data to generate 365 day data

library(tidyverse)
library(data.table) 
library(lubridate) 

#Load data
#Load Data: saved in two chunks currently
results1 <- readRDS("EMEP4uk_results_part1.rds")
results2 <- readRDS("EMEP4uk_results_part2.rds")

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
time_pd <- days(183)
app_dates$appt <- as.Date(app_dates$appt)
app_dates1 <- app_dates %>%
  mutate(start_date = appt - time_pd,
         id = as.character(id))
head(app_dates1)
sum(is.na(app_dates1))
app_dates1 <- app_dates1 %>%
  rename(ID = id) %>%
  mutate(ID = as.integer(ID))
  dim(app_dates1) 

#add in dates to results df
full_results <- left_join(full_results, app_dates1, by = "ID")

#check methylation IDs
meth_ids <- fread("2024-10-31_meth_ids.csv")
#Load list of IDs to keep
gs_ids <- read.table("GS_ids_20241007.txt") 
full_resultsids <- unique(full_results$ID)
length(full_resultsids) 
sum(is.na(full_results))


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
#filter to 183 days before the appointment (not including appointment day)
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
filtered_results[ID == "33", .N, by = pollutant] #183 entries for each pollutant 

# Count NAs for each column
na_counts <- filtered_results[, lapply(.SD, function(x) sum(is.na(x)))]
na_counts_full <-  full_results[, lapply(.SD, function(x) sum(is.na(x)))]

#*** compare don't overwrite
#Save out full183day results
saveRDS(filtered_results, file = "daily_183_pollution.RDS")  

#_________________________________________________________________________________________
#Create summary data

#add year column
filtered_results <- filtered_results %>%
 mutate(Year = year(measurement_day))

#Create 365day average data
annual_df <- filtered_results[, .(
  days_183_avg = mean(Measurement, na.rm = TRUE),
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
annual_df3 <- annual_df2 %>% drop_na()


#Save out scaled
annual_df4 <- annual_df3 %>%
  dplyr::select(ID, Pollutant, days_183_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_183_avg
  )
scaled_pollution <- annual_df4  
scaled_pollution[2:9] <- scale(annual_df4[2:9])   

#*** save out
fwrite(annual_df3, file = "6month_avg_exp_foranalysis_ntr.csv")
fwrite(scaled_pollution, file = "/6month_avg_exp_foranalysis_scaled.csv")



