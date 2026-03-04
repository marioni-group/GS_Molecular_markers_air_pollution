#Generating average pollutant data for participants from 2005 - 2011
#Load libraries
library(stringr)
library(purrr)
library(tictoc)
library(data.table)
library(tidyverse)
library(lubridate)
library(egg)
library(here)

#Load Data: saved in two chunks currently
results1 <- readRDS("/Pollution/EMEP4uk_gs/data_2005_2011/EMEP4uk_results_part1.rds")
results2 <- readRDS("/Pollution/EMEP4uk_gs/data_2005_2011/EMEP4uk_results_part2.rds")

#create dataframe of results
df <- dplyr::bind_rows(results1)
df2 <- dplyr::bind_rows(results2)

#combine together
full_results <- rbind(df, df2)
gc()
#_________________________________________________________________________________________________________
#3. Distribution of values by year & pollutant

#Load complete data (with IDs) - see above

#create yearly average data
full_results <- full_results %>%
  mutate(Date = ymd(Date))
full_results <- full_results %>% #extract Year and months from daily data
  mutate(Year = year(Date),
         Month = month(Date))

yearly_df <- full_results %>%
  group_by(ID, pollutant, Year) %>%  # Group by ID, pollutant, and year
  summarise(yearly_avg = mean(Measurement, na.rm = TRUE),
            standard_dev = sd(Measurement, na.rm = TRUE)) %>%  # Calculate yearly average
  ungroup()  # Ungroup to return a regular data frame
head(yearly_df)

fwrite(yearly_df, file = here("Pollution/EMEP4uk_gs/yearly_avg_exp.csv"))

#Create 7 year average data
annual_df <- full_results %>%
  group_by(ID, pollutant) %>%  # Group by ID, pollutant
  summarise(avg_7yrs = mean(Measurement, na.rm = TRUE),
            standard_dev = sd(Measurement, na.rm = TRUE),
            year_range = list(unique(Year))) %>%  # 
  ungroup()  # Ungroup to return a regular data frame

sapply(annual_df, function(y) sum(length(which(is.na(y)))))
head(full_results)
sapply(full_results, function(y) sum(length(which(is.na(y)))))

#________________________________________________________________________________________________
#Clean up pollution exposure df
#Complete pollutant abbreviations
annual_df <- annual_df %>% rename(raster_entry = pollutant)
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
fwrite(annual_df2, file = here("Pollution/EMEP4uk_gs/sevenyr_foranalysis_notsc.csv"))

#Filter to IDs in 365yr data
length(unique(annual_df2$ID))
length(unique(oneyr$ID))
annual_df2 <- annual_df2 %>% dplyr::filter(ID %in% oneyr$ID)
length(unique(annual_df2$ID)) 


#Scale
annual_dft <- annual_df2 %>% dplyr::select(ID, year_range, avg_7yrs, Pollutant)
annual_dft <- annual_dft %>%
  tidyr::pivot_wider(names_from = Pollutant, values_from = avg_7yrs)
annual_dft$ID <- as.character(annual_dft$ID)
annual_dft2 <- annual_dft %>% drop_na()
dim(annual_dft2)

annual_dfts <- annual_dft2 %>%
  dplyr::select(-year_range) %>%
   modify_if(is.numeric, ~scale(., center = TRUE, scale = TRUE))
fwrite(annual_dfts, file = "Pollution/EMEP4uk_gs/sevenyr_foranalysis_sc.csv")


#Create summary table for 7-year data
#Create a summary table
library(arsenal)
  
annual_df2 <- annual_df2 %>%
  mutate(year_range2 = c("2005:2011"))  
unique(annual_df2$year_range)


control <- tableby.control(,    # Use ANOVA for numeric tests
  total = TRUE,              # Include a total column
  numerics = list(
    mean = FALSE,
    Range = FALSE,
    median = TRUE,           # Show median
    q1q3 = TRUE)) 

tab1 <- summary(tableby(Pollutant ~ avg_7yrs, data = annual_df2, control = control), text = TRUE)
fwrite(tab1, file = "Pollution/EMEP4uk_gs/summary_table_7yrdata.csv")

