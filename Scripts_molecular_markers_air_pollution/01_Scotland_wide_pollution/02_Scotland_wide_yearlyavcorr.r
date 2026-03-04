#Exploratory pollution analysis: Calculate correlations between years

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


#______________________________________________________________________________
#Correlations using Terra -> plots

library(terra)
library(dplyr)
library(stringr)

pollutants <- c("SURF_ug_NO2", "SURF_ug_NO", "SURF_ug_PM25_rh50", "SURF_ppb_O3", "SURF_ug_PM10_rh50", "SURF_ug_SO2", "SURF_ug_NO3_C", "SURF_ug_NO3_F")
years <- as.character(c(2005:2011))
file_list <-  list.files(path = "Pollution/EMEP4UK_data", pattern = "*.nc", full.names = TRUE)

yearly_correlation <- function(years, pollutant, file_list) {
  message("Processing ", pollutant)

  scotland_map <- ne_states(geounit = "scotland")
  scotland_vector <- vect(scotland_map)

  # Extract files for each year
  files <- lapply(years, function(year) {
    file <- file_list[str_detect(file_list, year)]
    if (length(file) == 0) stop(paste("No files found for year:", year))
    file
  })

  # Create and project all raster stacks
  stacks <- lapply(files, function(file) terra::rast(file, subds = pollutant))
  scotland_vector <- project(scotland_vector, stacks[[1]])

  # Crop and calculate means
  mean_rasters <- lapply(stacks, function(stack) {
    cropped <- crop(stack, scotland_vector, snap = "out", mask = TRUE)
    app(cropped, mean, na.rm = TRUE)
  })

  # Combine all mean rasters
  combined <- do.call(c, mean_rasters)
  combined <- na.omit(combined)

  # Calculate the correlation matrix
  correlation <- layerCor(combined, cor, use = "complete.obs")
  rownames(correlation) <- years
  colnames(correlation) <- years

  return(correlation)
}


tic()
results <- map(pollutants, ~ yearly_correlation(years, .x, file_list))
toc()

names(results) <- pollutants
head(results)

saveRDS(results,file = "/Pollution/spatial_corr.rds")

#Make heatmaps

#test
test_map <- as.matrix(results[[1]])
Heatmap(test_map,
  row_order = rownames(test_map),
      column_order = colnames(test_map))


make_heatmap2 <- function(x) {
  subset <- results[[x]]
  ht <- Heatmap(as.matrix(subset),
            name = x,
      row_order = rownames(subset),
      column_order = colnames(subset),
  col = colorRamp2(c(global_min, (global_min + global_max) / 2, global_max), c("blue", "white", "red")),
  column_title = "",
  row_title = "",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),)
}

# Calculate overall min and max values across all results
all_values <- unlist(lapply(results, function(res) as.matrix(res)))
global_min <- min(all_values, na.rm = TRUE)
global_max <- max(all_values, na.rm = TRUE)



heatmap_list <- purrr::map(pollutants, make_heatmap2) %>%
  purrr::set_names(pollutants)
convert_togrob <- function(x) {
  grid.grabExpr(draw(x))
}  
heatmap_grobs <- purrr::map(heatmap_list, convert_togrob) 

library(gridExtra)

pdf(file = "Pollution/spatial_corr_htmaps.pdf")
grid.arrange(grobs = heatmap_grobs, widths = c(2,2), ncol = 2, labels = pollutants)
dev.off()

#-______________________________________________________________________________
#Calculating pairwise correlations using cor.test -> to generate full stats for suppl table

results_df <- fread("Pollution/EMEP4UK_yearlysummary.csv")

results_df <- results_df %>% dplyr::rename(raster_entry = "Pollutant")
results_df2 <- results_df %>%
  mutate(Pollutant = ifelse(raster_entry == "SURF_ppb_O3", "O3", 
                      ifelse(raster_entry == "SURF_ug_NO2", "NO2", 
                        ifelse(raster_entry == "SURF_ug_PM10_rh50", "PM10",
                          ifelse(raster_entry == "SURF_ug_PM25_rh50", "PM25",
                            ifelse(raster_entry == "SURF_ug_SO2", "SO2", 
                              ifelse(raster_entry == "SURF_ug_NO", "NO", 
                      ifelse(raster_entry == "SURF_ug_NO3_C", "NO3_C", 
                        ifelse(raster_entry == "SURF_ug_NO3_F", "NO3_F", raster_entry)))))))))
pollutant_list <- unique(results_df2$Pollutant)



results4 <- results_df2 %>%
  pivot_wider(
    names_from = Year,
    values_from = mean
  )
head(results4)  

# Compare pairwaise years 
years <- as.character(unique(results_df$Year))
year_pairs <- tidyr::expand_grid(years, years)
year_pairs <- year_pairs %>% dplyr::rename(year1 = "years...1", year2 = "years...2")
pair_labels <- year_pairs %>% 
   mutate(labels = paste(year1, "-", year2)) 

year_pair_corr <- function(b) {
  subset <- results4 %>% dplyr::filter(Pollutant == b)
test <- 
  purrr::map2(year_pairs$year1, year_pairs$year2, ~ {
  tidy(cor.test(subset[[.x]], subset[[.y]]))
})
corr_results <- bind_rows(test)
corr_results$years <- pair_labels$labels
corr_results$Pollutant <- b
return(corr_results)
}

df <- purrr::map(pollutant_list, year_pair_corr)
names(df) <- pollutant_list
by_yearcorr_df <- bind_rows(df, .id = "Pollutant")
by_yearcorr_df
fwrite(by_yearcorr_df, file = "Pollution/EMEP4UK_yearlycorrelations.csv")


#Plot to compare to previously generated heatmap
#Make a matrix
year_year_matrix <- function(d) {
  subset <- results4 %>% dplyr::filter(Pollutant == d)
  corr_matrix <-
    purrr::map2(year_pairs$year1, year_pairs$year2, ~
      cor.test(subset[[.x]], subset[[.y]])$estimate) %>%
   matrix(nrow = length(years), 
                             ncol = length(years),
                             dimnames = list(years, years))
}


  matrix_list <- purrr::map(pollutant_list, year_year_matrix) %>%
  purrr::set_names(pollutant_list)
head(matrix_list)

x <- matrix_list[[1]]
name <- name_list[[1]]
x <- as.matrix(x)
x[] <- as.numeric(x[])

make_heatmap3 <- function(x, name) {
  ht <- Heatmap(x,
            name = name,
      row_order = rownames(x),
      column_order = colnames(x),
  col = circlize::colorRamp2(c(-1, 0, 1), c("green", "grey", "red")),
  column_title = "",
  row_title = "",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8))
  return(ht)
}

# Calculate overall min and max values across all results
all_values <- unlist(lapply(matrix_list, function(res) as.matrix(res)))
global_min <- min(all_values, na.rm = TRUE)
global_max <- max(all_values, na.rm = TRUE)

name_list <- names(matrix_list)

 heatmap_list <- purrr::map(matrix_list, names(matrix_list),  make_heatmap3) 
#   purrr::set_names(name_list)
convert_togrob <- function(x) {
  grid.grabExpr(draw(x))
}  
heatmap_grobs <- purrr::map(heatmap_list, convert_togrob) 

library(gridExtra)
grid.arrange(grobs = heatmap_grobs, widths = c(2,2), ncol = 2, labels = pollutants) #This is the same as the spatial correlation


#____________________________________________________________________________________-
#Aggregate year correlations

#Create aggregate variables
results5 <- results4 %>%
  rename_with(~ paste0("year_", .), .cols = 5:11)
head(results5)

aggregate_year_summary <- results5 %>%
  dplyr::group_by(Pollutant) %>%
  summarise(year_2005 = mean(year_2005),
            year_2005_6 = mean(year_2005 + year_2006),
            year_2005_7 = mean(year_2005 + year_2006 + year_2007),
            year_2005_8 = mean(year_2005 + year_2006 + year_2007 + year_2008),
        year_2005_9 = mean(year_2005 + year_2006 + year_2007 + year_2008 + year_2009),
        year_2005_10 = mean(year_2005 + year_2006 + year_2007 + year_2008 + year_2009 + year_2010),
        year_2005_11 = mean(year_2005 + year_2006 + year_2007 + year_2008 + year_2009 + year_2010 + year_2011))


aggregate_years <- results5 %>%
  dplyr::group_by(Pollutant) %>%
  rowwise() %>%
  mutate(year_2005_6 = mean(c(year_2005, year_2006)),
        year_2005_7 = mean(c(year_2005, year_2006, year_2007)),
        year_2005_8 = mean(c(year_2005, year_2006, year_2007, year_2008)),
        year_2005_9 = mean(c(year_2005, year_2006, year_2007, year_2008, year_2009)),
        year_2005_10 = mean(c(year_2005, year_2006, year_2007, year_2008, year_2009, year_2010)),
        year_2005_11 = mean(c(year_2005, year_2006, year_2007, year_2008, year_2009, year_2010, year_2011)))



results6 <- aggregate_years %>% dplyr::select(!c("year_2006", "year_2007", "year_2008","year_2009","year_2010","year_2011"))

aggr_years <- colnames(results6)[5:11]
aggr_pairs <- tidyr::expand_grid(aggr_years, aggr_years)
aggr_pairs <- aggr_pairs %>% dplyr::rename(year1 = "aggr_years...1", year2 = "aggr_years...2")
pair_labels <- aggr_pairs %>% 
   mutate(labels = paste(year1, "-", year2)) 

agg_pair_corr <- function(b) {
  subset <- results6 %>% dplyr::filter(Pollutant == b)
test <- 
  purrr::map2(aggr_pairs$year1, aggr_pairs$year2, ~ {
  tidy(cor.test(subset[[.x]], subset[[.y]], use = "complete.obs"))
})
corr_results <- bind_rows(test)
corr_results$years <- pair_labels$labels
corr_results$Pollutant <- b
return(corr_results)
}

df <- purrr::map(pollutant_list, agg_pair_corr)
by_aggrcorr_df <- bind_rows(df, .id = "Pollutant")
fwrite(by_aggrcorr_df, file = "Pollution/EMEP4UK_aggr_yearlycorrelations.csv")


