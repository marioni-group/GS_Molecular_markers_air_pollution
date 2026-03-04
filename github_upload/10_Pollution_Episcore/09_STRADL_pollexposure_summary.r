#Process pollution exposure 

# sh-4.4$ eval "$(micromamba shell hook --shell bash)"
# sh-4.4$ micromamba activate R

library(tidyverse)
library(data.table) 
library(lubridate)

#Load data
results <- readRDS("EMEP4uk_results_STRADL.rds")

#create dataframe of results
df <- dplyr::bind_rows(results)
(df)
df <- df %>% mutate(Date = lubridate::as_date(Date))
df <- df %>% 
  mutate(measurement_year = lubridate::year(Date))

#Create 365day average data
df <- setDT(df)
annual_df <- df[, .(
  days_365_avg = mean(Measurement, na.rm = TRUE),
  standard_dev = sd(Measurement, na.rm = TRUE),
  year_range = list(unique(measurement_year))
), by = .(ID, pollutant)]

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
unique(annual_df2$year_range)
fwrite(annual_df2, file = "stradl_365_av_notrs.csv")

#scale
sum(is.na(annual_df2))
annual_df3 <- annual_df2 %>%
  dplyr::select(ID, Pollutant, days_365_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )

scaled_pollution <- annual_df3  
scaled_pollution[2:9] <- scale(annual_df3[2:9])   
str(scaled_pollution)
fwrite(scaled_pollution, file = "stradl_365_av_sc.csv")
