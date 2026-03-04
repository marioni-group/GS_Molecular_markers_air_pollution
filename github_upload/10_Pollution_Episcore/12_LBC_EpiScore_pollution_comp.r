#LBC_EpiScore_pollution_comp

#Compare EpiScore and pollution exposure

library(tidyverse)
library(data.table)

lbc <- fread("lbc_wv1_episcores.csv")
pollution <- fread("lbc_365day_av_101025.csv")

#tidy up pollution data
pollution <- pollution %>%
  mutate(Pollutant = ifelse(raster_entry == "SURF_ug_NO", "NO", 
                      ifelse(raster_entry == "SURF_ug_NO3_C", "NO3_C", 
                        ifelse(raster_entry == "SURF_ug_NO3_F", "NO3_F",
                          ifelse(raster_entry == "SURF_ug_PM25_rh50", "PM25",
                            ifelse(raster_entry == "SURF_ppb_O3", "O3", 
                      ifelse(raster_entry == "SURF_ug_NO2", "NO2", 
                        ifelse(raster_entry == "SURF_ug_PM10_rh50", "PM10",
                            ifelse(raster_entry == "SURF_ug_SO2", "SO2", Pollutant)))))))))

pollution <- pollution %>%
  dplyr::select(LBC36no, Pollutant, days_365_avg) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = days_365_avg
  )

scaled_pollution <- pollution
scaled_pollution[2:9] <- scale(pollution[2:9])   
str(scaled_pollution)
fwrite(scaled_pollution, file = "lbcwv1_365_av_sc.csv")



lbc <- dcast(melt(lbc, id.vars = "pollutant"), variable ~ pollutant)


episcore <- setDT(lbc)
episcore$variable <- as.character(episcore$variable)

#Match IDs and orders
wv1 <- fread("wv1_target36.csv")

#Go from lbcno to sample basename
wv1 <- wv1 %>% select(ID_raw, Basename)
wv1 <- wv1 %>% rename(LBC36no = "ID_raw")
poll <- left_join(scaled_pollution, wv1, by = "LBC36no")

#filter to match entries
poll <- poll %>% filter(Basename %in% episcore$variable)


#scale episcore
str(episcore)
episcore <- episcore %>% filter(variable %in% poll$Basename)
episcore <- episcore %>% mutate(across(where(is.numeric), scale))

#match orders of pollution & episcore
episcore <- episcore[match(poll$Basename, episcore$variable),]
identical(episcore$variable, poll$Basename) 

#Prep for plotting/regression
episcore <- episcore %>% mutate(type = "episcore") %>% rename(Basename = variable)
episcore <- episcore %>% pivot_longer(.,
    !c(Basename, type),
    names_to = "pollutant",
    values_to = "value")

measured <- poll %>% mutate(type = "modelled") %>% select(-LBC36no)
measured <- measured %>% pivot_longer(.,
     !c(Basename, type),
     names_to = "pollutant",
     values_to = "value")

#combine
plot_df <- rbind(episcore, measured)
plot_df1 <- plot_df %>% pivot_wider(.,
                                    names_from = "type",
                                   values_from = "value")
head(plot_df1)                                   
sum(is.na(plot_df1))
fwrite(plot_df1, file = "lbc36wv1_modelled_episcore.csv")


#Function calculating r
cor_predictr <-   function(x){
  return(
          broom::tidy(cor.test(x$episcore, x$modelled, method = "pearson"))) 
}

#on 693 individuals
correlation_results <- plot_df1 %>%
  filter(!pollutant == "logsmok") %>%
  group_by(pollutant) %>%
  do(cor_predictr(.))
correlation_results <- correlation_results %>% arrange(desc(estimate))
correlation_results

fwrite(correlation_results, file = "cor_lbcwv1.csv")
res <- fread("cor_lbcwv1.csv")

#_____________________________________________________
#Plot
#_____________________________________________________

level_order <- correlation_results$pollutant


#PLOT.
p1 <- correlation_results %>%
  ggplot(aes(x = factor(pollutant, level = level_order), y = estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  theme_bw() +
  labs(x = "Pollutant EpiScore", y = "Cor. with pollution exposure") +
  theme(axis.title.x = element_text(margin = margin(t = 30), size = 10),
        axis.title.y = element_text(margin = margin(r = 20), size = 10),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 8)) +
  coord_fixed(ratio = 15)
ggsave("lbc_episcore_poll_corr.pdf")


