#LBC Geocoding

library(tidyverse)
library(data.table)
library(haven)
library(terra)
library(PostcodesioR)
library(purrr)
library(lubridate)
library(tictoc)

#Load data
lbcgeo <- read_sav("lbc_dat_w1.sav")
head(lbcgeo)

#__________________________________________________________
#1.Create latitude and longitude for list of postcodes
#__________________________________________________________

#Generate latitude and longitude for postcodes
#Extract latitude and longitude using UK postcode lookup
extract_latlong <- function (x) {
postcode_info <- x %>% purrr::map(postcode_lookup) #extracts data for each postcode
postcode_latlong <- purrr::map(postcode_info, ~dplyr::select(., postcode, longitude, latitude)) 
latlong_df <- bind_rows(postcode_latlong) 
return(latlong_df)
}

#Extract latitude and longitude for Scottish postcode lookup only
extract_latlong_scottish <- function (x) {
postcode_info <- x %>% purrr::map(scottish_postcode_lookup) #extracts data for each postcode
latlong_df <- bind_rows(postcode_info) 
return(latlong_df)
}




postcodes <- lbcgeo$postcode_w1
latlong_df <- extract_latlong(postcodes) #3 postcodes not available/found

#list of postcodes no longer current
not_current <- latlong_df %>% filter(is.na(longitude)) %>% select(postcode)
not_current <- not_current %>% mutate(postcode = gsub("%20", " ", postcode))
not_current <- not_current$postcode

head(latlong_df)
latlong_df2 <- extract_latlong_scottish(postcodes) #if use scottish look-up don't get latitude/longitude, highlights 3 postcodes which are not scottish
latlong_df2 <- latlong_df2 %>% mutate(postcode = gsub("%20", " ", postcode))
not_scottish <- latlong_df2 %>% filter(is.na(scottish_parliamentary_constituency)) %>% select(postcode)
not_scottish <- not_scottish %>% filter(!postcode %in% not_current)
not_scottish <- not_scottish$postcode

#create df of GeoLocID, latitude/longitude, flag for non-scottish
latlong_df2 <- latlong_df2 %>% 
  mutate(flag = ifelse(is.na(scottish_parliamentary_constituency) & postcode %in% not_scottish, "3",  #postcodes removed
                  ifelse(is.na(scottish_parliamentary_constituency) & postcode %in% not_current, "2", "1")))
x <- latlong_df2 %>% select(postcode, flag)
x <- distinct(x)
x <- x %>% mutate(postcode = gsub("%20", " ", postcode))
x %>%
  group_by_all() %>%
  filter(n()>1) %>%
  ungroup()

#remove duplicate postcodes
latlong_df <- distinct(latlong_df)
latlong_df <- latlong_df %>% mutate(postcode = gsub("%20", " ", postcode))

latlong_df3 <- left_join(latlong_df, x, by = "postcode") 
latlong_df3 %>% #check no duplicates
  group_by_all() %>%
  filter(n()>1) %>%
  ungroup()

               
lbcgeo <- lbcgeo %>% rename(postcode = postcode_w1) %>% select(GeoLocID, postcode)
df <- left_join(lbcgeo, latlong_df3, by = "postcode")
dim(df)
head(df)

#Prep for saving out
df <- df %>% select(-postcode)

df %>%
  group_by_all() %>%
  filter(n()>1) %>%
  ungroup()


fwrite(df, file = "lat_longdf_081025.csv")




#Reload data
lbcgeo <- read_sav("lbc_dat_w1.sav")

#Merge to date info for pollution data extraction
lbcgeo <- lbcgeo %>% rename(postcode = postcode_w1)
lbcgeo2 <- left_join(lbcgeo, latlong_df, by = "postcode")
head(lbcgeo2)
head(latlong_df)

#for this project - remove flagged postcodes
lbcgeo2 <- lbcgeo2 %>% filter(!postcode %in% not_current & !postcode %in% not_scottish)


#__________________________________________________________
#2. Generate 365 day time window pre baseline appointment
#__________________________________________________________

#Create an interval of 365 days
time_pd <- days(365)
lbcgeo2$datetested_w1 <- as_date(lbcgeo2$datetested_w1)

lbcgeo2 <- lbcgeo2 %>%
  mutate(start_date = datetested_w1 - time_pd,
         GeoLocID = as.character(GeoLocID))


#__________________________________________________________
# Extract data - by pollutant:
    #Create list of raster files 
    #Create list of locations for extraction & convert to correct projection
    #Turn lon/lat into points SpatVector
    #Create list of pollutants to extract data for
    #Extract pollutant data from raster
    #Extract data for each individual
#__________________________________________________________

#load file list - this folder contains raster files for relevant years of EMEP4uk model
file_list <-  list.files(path = "/EMEP4uk", pattern = "*.nc", full.names = TRUE)
#update to correct years (2003 - 2007)
file_list <- file_list[4:8]

#Create locations for data extraction
#set latitude and longitudes for the individuals (see prev script for latlong_df2 generation)
lon <- lbcgeo2$longitude
lat <- lbcgeo2$latitude
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
latlong_ids <- lbcgeo2$GeoLocID
start_dates <- lbcgeo2$start_date
appt_dates <- lbcgeo2$datetested_w1


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

    #add units for pollutant
    pollutant_units <- unique(units(r))

        poll_data <- purrr::map(latlong_ids, extract_polldata)
    poll_df <- bind_rows(poll_data)
    EMEP4uk_results[[pollutant]] <- poll_df
}
toc()

saveRDS(EMEP4uk_results, file = "/LBC_emep4uk_081025.rds")
full_results <- bind_rows(EMEP4uk_results, .id = "pollutant")
head(full_results)
fwrite(full_results, file = "/LBC_emep4uk_081025.csv")


res <- fread("/LBC_emep4uk_081025.csv")
head(res)
length(unique(res$ID))
length(unique(res$pollutant))


check <- res %>% count(ID, pollutant)
head(check)
check %>% filter(ID == "1")
check %>% filter(!n == 365) #all categories have 365 entries

