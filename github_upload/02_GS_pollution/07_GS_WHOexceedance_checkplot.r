#WHO limit exceedances
#Cut-offs from: https://www.who.int/news-room/feature-stories/detail/what-are-the-who-air-quality-guidelines

library(data.table)
library(tidyverse)

#set up WHO cut-offs
#Annual and 24hr cut-offs

WHO_limits <- data.frame(
  Pollutant = c("PM25", "PM25", "PM10", "PM10", "O3","O3", "NO2","NO2", "SO2", "CO"),
  Averaging_time = c("Annual", "24hr", "Annual", "24hr", "Peak_season", "8hr", "Annual", "24hr", "24hr", "24hr"),
  AQG_2021 = c(5,15,15,45,60,100, 10,25,40,4)
)
fwrite(WHO_limits, file = here("who_limits.csv"))


#Load daily pollutant-exposure for 365-day exposure
daily_avg <- readRDS("daily_365_pollution.RDS")

#Filter to pollutants of interest
#add pollutant abbrev. for easier data manipulation
daily_avg <- daily_avg %>%
  mutate(Pollutant = ifelse(
    pollutant == "SURF_ug_NO2", "NO2", ifelse(
      pollutant == "SURF_ug_PM10_rh50", "PM10", ifelse(
        pollutant == "SURF_ppb_O3", "O3", ifelse(
          pollutant == "SURF_ug_SO2", "SO2", ifelse(
            pollutant == "SURF_ug_PM25_rh50", "PM25", "not_available")
          )
        )
      )
    )
  )

#Rename for clarity
daily_avg <- daily_avg %>%
  rename(raster_entry = pollutant)

#Subset to pollutants of interest
who_poll <- unique(WHO_limits$Pollutant)

# Convert for efficiency
setDT(daily_avg)
library(foreach)
library(doParallel)

# Set up cluster
cl <- makeCluster(12)
registerDoParallel(cl)

# Parallel processing with foreach
daily_avg2 <- foreach(split_data = iter(split(daily_avg, by = "ID")), .packages = "data.table") %dopar% {
  split_data[Pollutant %in% who_poll]
}

# Stop cluster
stopCluster(cl) 

# Combine the results
daily_avg2 <- rbindlist(daily_avg2)

#Calculate exceedances by day 
WHO_limits_24hr <- WHO_limits %>% dplyr::filter(Averaging_time == "24hr")
daily_avg2 <- left_join(daily_avg2, WHO_limits_24hr, by = "Pollutant", relationship = "many-to-many")  
daily_avg2 <- daily_avg2 %>%
  # dplyr::filter(Averaging_time == "24hr") %>%
  mutate(WHO_limit_exceedance = ifelse(
    Measurement > AQG_2021, "yes", "no"    
  ))

#check IDs (should be 22,071)

#Summarise
daily_exceedances <- daily_avg2 %>% 
  group_by(ID, Pollutant) %>%
  summarise(Daily_exceedances = sum(Measurement > AQG_2021)) %>%
  ungroup()
fwrite(daily_exceedances, file = "daily_exceedances_who.csv")

#Plot
library(RColorBrewer)
col <- RColorBrewer::brewer.pal(9, 'Paired')
col <- col[1:2]

#*****check plot****** - saved in datastore = correct, saved on cluster = not correct
pdf("daily_exceedances_who.pdf")
daily_exceedances %>%
  dplyr::filter(!Pollutant == "O3") %>%
  mutate(above = ifelse(Daily_exceedances > 4, "yes", "no")) %>%
  ggplot(aes(x = Daily_exceedances, fill = above)) +
  geom_histogram(bins = 100) +
  scale_fill_manual(values = col) +
  xlab("Number of days exceeding 24hr limit") +
  ylab("Participant count") +
  labs(fill = "> 4 exceedance days") +
  theme_bw() +
  facet_wrap(~Pollutant, scales = "free_y") 
dev.off()  


#Calculate overall exceedances - according to WHO guidelines 3-4 exceedance days are allowable
#filter to IDs in analysis 
setDT(daily_exceedances)
daily_exceedances2 <- daily_exceedances[ID %in% annual_avg$ID]

overall_exceedances <- daily_exceedances2 %>%
  dplyr::filter(!Pollutant == "O3") %>%
  mutate(above = ifelse(Daily_exceedances > 4, TRUE, FALSE)) %>% 
  group_by(Pollutant) %>%
  summarise(overall_exceedances = sum(above, na.rm = TRUE)) %>%
  ungroup()
fwrite(overall_exceedances, file = "overall_exceedances.csv")

### annual average exposure_________________________________________________________

#Load average pollutants for 365-day exposure
annual_avg <- fread("annual_avg_exp_foranalysis_ntr.csv")

annual_avg2 <- left_join(annual_avg, WHO_limits, by = "Pollutant", relationship = "many-to-many")
# #For each pollutant - count number exceedances of WHO limit for annual exposure

#Using 365 average exposure data
annual_exceedances <- annual_avg2 %>%
  drop_na() %>%
  dplyr::filter(Averaging_time == "Annual") %>%
 mutate(WHO_limit_exceedance = ifelse(days_365_avg > AQG_2021, "yes", "no"))

 #save out
 fwrite(annual_exceedances, file = "365day_ex_wholimits.csv")

annual_sum_forplot <- annual_exceedances %>%
  group_by(Pollutant) %>%
  summarise(annual_exceedances = sum(days_365_avg > AQG_2021)) #SO2 does not have an annual limit
fwrite(annual_sum_forplot, file = "annual_exp_exceedances.csv")  

