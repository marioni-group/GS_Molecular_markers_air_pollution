#Script to extract pollution data for GS addresses

library(ncdf4)
library(raster)
library(lattice)
library(PostcodesioR)
library(RMariaDB)
library(dplyr)
library(sf)
library(sgo)
library(purrr)
library(lubridate)
library(tidyverse)
library(data.table)
library(terra)
library(tictoc)


#__________________________________________________________
#1.Create latitude and longitude for list of postcodes
#__________________________________________________________

#Using dummy postcodes for now
postcodes <- c("") #STRADL postcodes

#Generate latitude and longitudes for postcodes
#Extract information for each postcode in list using postcodioR package
extract_latlong <- function (x) {
postcode_info <- x %>% purrr::map(postcode_lookup) #extracts data for each postcode
postcode_latlong <- purrr::map(postcode_info, ~dplyr::select(., postcode, longitude, latitude)) 
latlong_df <- bind_rows(postcode_latlong) 
return(latlong_df)
}

latlong_df <- extract_latlong(postcodes)
latlong_df



latlong_df <- left_join(latlong_df, ids, by = "postcode")

#__________________________________________________________
#2 & 3. Find baseline appointment date for linked ID & generate 365 day time window
#__________________________________________________________

#use STRADL appointment dates as end date for pollution capture e.g. 
app_dates <- fread("") #add in STRADL appointment dates
app_dates$appt <- gsub(" .*", "", app_dates$appt) #remove time from appt as don't require this
head(app_dates)

#Create an interval of 365 days
time_pd <- days(365)
app_dates$appt <- as.Date(app_dates$appt)
app_dates1 <- app_dates %>%
  mutate(start_date = appt - time_pd,
         id = as.character(id))

#join with latlong_df
latlong_df2 <- left_join(latlong_df, app_dates1, by = "id")


#__________________________________________________________
# Extract data - by pollutant:
    #Create list of raster files 
    #Create list of locations for extraction & convert to correct projection
    #Turn lon/lat into points SpatVector
    #Create list of pollutants to extract data for
    #Extract pollutant data from raster
    #Extract data for each individual
#__________________________________________________________
#load file list
file_list <-  list.files(path = "EMEP4UK_data", pattern = "*.nc", full.names = TRUE)
file_list #should display files from 2014 - 2017

#Create locations for data extraction
#set latitude and longitudes for the individuals (see prev script for latlong_df2 generation)
lon <- latlong_df2$longitude
lat <- latlong_df2$latitude
pts1 <- cbind(lon,lat)
pts_sv <- terra::vect(pts1, type = "points", crs = "EPSG:4326") #lat and lon are of this crs

#load one raster to get projection
r <- terra::rast(file_list[1], subds = "SURF_ug_PM25_rh50") 

#reproject point data to match crs of r file
pts_sv2 <- terra::project(pts_sv, r, partial = FALSE)

#select pollutants for extraction
pollutant_list <- c("SURF_ug_PM25_rh50", "SURF_ug_NO2", "SURF_ug_NO", "SURF_ppb_O3", "SURF_ug_PM10_rh50", "SURF_ug_SO2", "SURF_ug_NO3_C", "SURF_ug_NO3_F") 



#Function to extract data for each ID 
extract_polldata <- function(i) {
            start <- start_dates[latlong_ids == i]
            end <- appt_dates[latlong_ids == i]

         # Subset columns based on the date range
        subset_cols <- layer_dates[layer_dates >= start & layer_dates < end]
        subset_cols <- as.character(subset_cols)

        poll <- extractedData %>% 
            dplyr::filter(ID == i) %>%
            dplyr::select(all_of(c("ID",subset_cols))) %>%
        pivot_longer(
                !ID,
                names_to = "Date",
                values_to = "Measurement") %>%
            mutate(
            pollutant = pollutant,
            units = pollutant_units)
        return(poll)
        }


EMEP4uk_results <- list() #initialise list

#Pre-extract latlong_df2 columns required
latlong_ids <- latlong_df2$id
start_dates <- latlong_df2$start_date
appt_dates <- latlong_df2$appt

tic()
for(pollutant in pollutant_list) {
    print(paste0("Working on pollutant", pollutant))
    r <- terra::rast(file_list, subds = pollutant) #creates raster of all years, for selected pollutant

    #extract data for longitude & latitude points
    extractedData <- terra::extract(r, pts_sv2, method = "bilinear")

    # Get the dates from the raster stack
    layer_dates <- as.Date(time(r))

    # Replace the column names with the corresponding dates
    colnames(extractedData)[-1] <- as.character(layer_dates)

    #Replace ID column with IDs from input df
    extractedData$ID <- latlong_ids

    #For each individual, subset to 365 days before their baseline appointment

    pollutant_units <- unique(units(r))

        poll_data <- purrr::map(latlong_ids, extract_polldata)
    poll_df <- bind_rows(poll_data)
    EMEP4uk_results[[pollutant]] <- poll_df
}
toc()
