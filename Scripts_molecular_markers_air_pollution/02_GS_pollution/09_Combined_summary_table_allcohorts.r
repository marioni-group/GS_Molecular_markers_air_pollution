#Combined summary table for assigned pollution data, for each time window in GS & LBC

#include: mean, median, IQR

#Function to summarise pollution data
stat_poll <- function (x) {
df <- data.frame( 
 Mean = mean(x),
 Median = median(x),
 Stdev = sd(x),
 q1 = unname(quantile(x, 0.25)),
 q3 = unname(quantile(x, 0.75)),
 iqr = IQR(x),
 Min = min(x),
 Max = max(x),
 row.names = NULL)
 return(df)
}



#GS 1year
oneyr <- fread("annual_avg_exp_foranalysis_ntr.csv")
oneyr <- oneyr %>% select(-c(year_range, standard_dev, raster_entry))
oneyr <- oneyr %>% pivot_wider(
    id_cols = NULL,
    names_from = Pollutant,
    values_from = days_365_avg
)


summary_dfmulti <- oneyr %>% 
  drop_na() %>%
  dplyr::select(-ID) %>%
  map(., stat_poll)

#Summary stats
summary_dfmulti <- bind_rows(summary_dfmulti, .id = "Pollutant")
summary_dfmulti1 <- summary_dfmulti %>% mutate(time_span = "oneyr")



#GS 6 month
sixmonth <- fread("6month_avg_exp_foranalysis_ntr.csv")

six_dfmulti <- sixmonth %>% 
  drop_na() %>%
  dplyr::select(-ID) %>%
  map(., stat_poll)

#Summary stats
six_dfmulti <- bind_rows(six_dfmulti, .id = "Pollutant")
six_dfmulti1 <- summary_dfmulti %>% mutate(time_span = "sixmonth")

#GS 7 year
annual_df2 <- fread("sevenyr_foranalysis_notsc.csv")
sevenyr <- annual_df2 %>% select(-c(year_range, standard_dev, raster_entry))
sevenyr <- sevenyr %>% pivot_wider(
    id_cols = NULL,
    names_from = Pollutant,
    values_from = avg_7yrs
)


seven_dfmulti <- sevenyr %>% 
  drop_na() %>%
  dplyr::select(-ID) %>%
  map(., stat_poll)

#Summary stats
seven_dfmulti <- bind_rows(seven_dfmulti, .id = "Pollutant")
seven_dfmulti1 <- seven_dfmulti %>% mutate(time_span = "sevenyr")


#combine all GS summaries
comb_summary <- rbind(summary_dfmulti1, six_dfmulti1, seven_dfmulti1)

#LBC 1 year
lbc_1year <- fread("lbc_365day_av_101025.csv")

lbc_1year <- lbc_1year %>% select(-c(units, standard_dev, raster_entry))
lbc_1year <- lbc_1year %>% pivot_wider(
    id_cols = NULL,
    names_from = Pollutant,
    values_from = days_365_avg
)

lbcsum <- lbc_1year %>% 
  drop_na() %>%
  dplyr::select(-LBC36no) %>%
  map(., stat_poll)

#Summary stats
lbcsum <- bind_rows(lbcsum, .id = "Pollutant")
lbcsum <- lbcsum %>% mutate(time_span = "lbc_oneyr")


#STRADL 1 year
STRADL_1year <- fread("stradl_365_av_notrs.csv")

stradl_1yr <- STRADL_1year %>% select(-c(year_range, standard_dev, raster_entry))
stradl_1yr <- stradl_1yr %>% pivot_wider(
    id_cols = NULL,
    names_from = Pollutant,
    values_from = days_365_avg
)

stradlsum <- stradl_1yr %>% 
  drop_na() %>%
  dplyr::select(-ID) %>%
  map(., stat_poll)

#Summary stats
stradlsum <- bind_rows(stradlsum, .id = "Pollutant")
stradlsum <- stradlsum %>% mutate(time_span = "stradl_oneyr")



allcomb <- rbind(comb_summary, lbcsum, stradlsum)
fwrite(allcomb, file = "poll_exps_combined_summary.csv")


#Calculate correlations for 1yr, 6month and 7year data in GS
#1yr - 7yr done in previous script

#6month - 1yr
head(sixmonth)
sixm <- sixmonth %>% pivot_longer(
    cols = !ID,
   names_to = "pollutant",
   values_to = "measurement" 
) %>%
dplyr::rename(sixmonth = measurement)


head(oneyr)
oney <- oneyr %>% pivot_longer(
    cols = !ID,
   names_to = "pollutant",
   values_to = "measurement" 
) %>%
dplyr::rename(oneyr = measurement)

identical(sixm$ID, oney$ID)
pollutant_list <- unique(oney$pollutant)
df <- full_join(sixm, oney, by = c("ID", "pollutant"))


corr_results <- list()
for (pollutant in pollutant_list) {
  x <- df %>% dplyr::filter(pollutant == pollutant)
  corr_results[[pollutant]] <- cor.test(x$oneyr, x$sixmonth)
}
corr_results


library(broom)
corr_results2 <- map(corr_results, tidy)
corr_results2
corr_results2 <- bind_rows(corr_results2, .id = "pollution")
head(corr_results2)
corr_results2 <- corr_results2 %>% mutate(inputs = "GS_oneyr_sixmonth")
oneyr_sixmnth <- corr_results2

#Then compare 6month - 7year
head(sevenyr)
seveny <- sevenyr %>% pivot_longer(
    cols = !ID,
   names_to = "pollutant",
   values_to = "measurement" 
) %>%
dplyr::rename(sevenyear = measurement)

#combine
dfc <- full_join(seveny, sixm, by = c("ID", "pollutant"))

corr_results <- list()
for (pollutant in pollutant_list) {
  x <- dfc %>% dplyr::filter(pollutant == pollutant)
  corr_results[[pollutant]] <- cor.test(x$sevenyear, x$sixmonth)
}
corr_results

corr_results2 <- map(corr_results, tidy)
corr_results2
corr_results2 <- bind_rows(corr_results2, .id = "pollution")
head(corr_results2)
corr_results2 <- corr_results2 %>% mutate(inputs = "GS_sevenyr_sixmonth")

#combine
corr_results <- rbind(oneyr_sixmnth, corr_results2)
fwrite(corr_results, file = "diff_timespan_corrs.csv")



#PLOTs

#. Compare 1year exposures for STRADL vs LBC vs GS
gs <- oneyr %>% mutate(cohort = "GS") %>% select(-ID)
stradl <- stradl_1yr %>% mutate(cohort = "STRADL") %>% select(-ID)
lbc <- lbc_1year %>% mutate(cohort = "LBC") %>% select(-LBC36no)

comb <- rbind(gs, stradl, lbc)
comb2 <- comb %>%
  pivot_longer(
    cols = !cohort,
    names_to = "pollutant",
    values_to = "exposure"
  )

library(paletteer)
combplot <- comb2 %>%
  ggplot(aes(x = cohort, y = exposure, fill = pollutant)) +
  geom_violin(trim = TRUE, draw_quantiles =  c(0.25, 0.5, 0.75), alpha = 0.8) + 
  scale_fill_paletteer_d("MetBrewer::Tam") +
  theme_minimal() +
  facet_wrap(~pollutant, scales = "free_y", nrow = 2)
ggsave("GS_STRADL_LBCcomp.pdf", height = 20, width = 25, units = "cm")
