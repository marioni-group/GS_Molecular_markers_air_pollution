#returned lbc data
#see scripts/pollution_analysis/LBC_geocoding/returned_lbc_data.r

library(tidyverse)
library(data.table)
library(haven)

full_results <- read_sav("LBC_emep4uk_081025_MERGED_CLEAN.sav")
full_results[1:5, 1:5]


#check data  are as expected
#____________________________________
setDT(full_results)

#should have 365 entries per pollutant per ID
check <- full_results[pollutant == "SURF_ug_PM10_rh50", .N, by = LBC36no]
na_counts <- full_results[, lapply(.SD, function(x) sum(is.na(x)))]
head(na_counts) #0
full_results[is.na(Measurement),] %>% select(LBC36no) %>% unique() #0
length(unique(full_results$LBC36no))

dim(full_results)
full_results2 <- distinct(full_results)
dim(full_results2) 
full_results2 <- drop_na(full_results2)
dim(full_results2) 


#Create 365day average data
annual_df <- full_results[, .(
  days_365_avg = mean(Measurement, na.rm = TRUE),
  standard_dev = sd(Measurement, na.rm = TRUE)
), by = .(LBC36no, pollutant)]

head(annual_df)
dim(annual_df)


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

#add in a units column
annual_df2 <- annual_df2 %>%
  mutate(units = ifelse(Pollutant == "O3", "ppb", "ug"))
head(annual_df2)
fwrite(annual_df2, file = "lbc_365day_av_101025.csv")

#___________________________________________
#Summarise data
#____________________________________________

#Table
library(arsenal)

control <- tableby.control(,    # Use ANOVA for numeric tests
  total = TRUE,              # Include a total column
  numerics = list(
    mean = FALSE,
    Range = FALSE,
    median = TRUE,           # Show median
    q1q3 = TRUE)) 
tab1 <- summary(tableby(Pollutant ~ days_365_avg, data = annual_df2, control = control), text = TRUE)
fwrite(as.data.frame(tab1), file = "lbc_365day_av_101025.csv")

tab <- fread("lbc_365day_av.csv")

#Boxplots
p1 <- annual_df2 %>%
 ggplot(aes(x = Pollutant, y = days_365_avg)) +
 geom_boxplot() + 
 theme_bw()
ggsave(p1, file = "summary_365avg_boxplots.pdf")

#correlations
df2 <- annual_df2 %>% dplyr::select(LBC36no, days_365_avg, Pollutant)
head(df2)

df2 |>
  dplyr::summarise(n = dplyr::n(), .by = c(Pollutant)) |>
  dplyr::filter(n > 1L) 

df3 <- df2 %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )

cor_poll <- cor(df3[2:9])
pdf("between_poll_htmp_365avg.pdf")
poll_heatmap <- ComplexHeatmap::Heatmap(cor_poll)
ComplexHeatmap::draw(poll_heatmap)
dev.off()

