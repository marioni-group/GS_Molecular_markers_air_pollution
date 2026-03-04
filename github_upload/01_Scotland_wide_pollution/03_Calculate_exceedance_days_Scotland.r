#03_Calculate_exceedance_days_Scotland
#Load libraries
library(terra)
library(raster)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(stringr)
library(purrr)
library(tictoc)
library(data.table)
library(tidyverse)
library(egg)
library(broom)

#Create file-list
file_list <- list.files(path = "Pollution/EMEP4uk", pattern = "*.nc", full.names = TRUE)

# Define list of pollutants and years and limits
pollutants <- c("SURF_ug_NO2", "SURF_ug_NO", "SURF_ug_PM25_rh50", "SURF_ppb_O3", "SURF_ug_PM10_rh50", "SURF_ug_SO2", "SURF_ug_NO3_C", "SURF_ug_NO3_F")
years <- as.character(c(2005:2011))
combinations <- tidyr::expand_grid(pollutants, years)
WHO_limits <- fread("Pollution/EMEP4uk_gs/who_limits.csv")

#load yearly summary
poll <- fread("Pollution/EMEP4uk/EMEP4UK_yearlysummary.csv") 

#Start with annual exceedances
#Function to generate number of exceedance days per raster per pollutant
yearly_exceedances <- function(pollutants, years) {
  message("Processing", pollutants, " for year ", years)

 #Extract correct file year
    file <- str_subset(file_list, as.character(years))
 #Load data for that year for selected pollutant   
    stack <- terra::rast(file, subds = pollutants)

 #crop to scotland
    scotland_map <- ne_states(geounit ="scotland")
    scotland_vector <- vect(scotland_map)
    crs(scotland_vector)
    crs(stack)
    scotland_vector <- project(scotland_vector, stack) #make sure in some coordinate system

 #create cropped raster
    cropped_stack <- crop(stack, scotland_vector, snap = "out", mask = TRUE)

 #create averaged raster
    output <- mean(cropped_stack, na.rm = TRUE)

  #Try converting to df for ggplot
    yearly_df <- as.data.frame(output, xy = TRUE)
    yearly_df <- yearly_df %>%
        mutate(Year = years,
                 Pollutant = pollutants)
  return(yearly_df)
}

# Apply a function to each combination
tic()
results <- pmap(combinations, yearly_means)
toc()

#Save: 
fwrite(results, file = "")