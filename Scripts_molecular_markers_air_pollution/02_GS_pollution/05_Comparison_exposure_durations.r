####Correlations
#Calculate between-pollutant correlations
df <- fread("annual_avg_exp_foranalysis_ntr.csv")
df2 <- df %>% dplyr::select(ID, days_365_avg, Pollutant)
head(df2)

df2 %>%
  dplyr::summarise(n = dplyr::n(), .by = c(Pollutant)) %>%
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

fwrite(cor_poll, file = "between_poll_corr_matrix.csv")


#______________________________________________________________________
#add in 6 month and 7 year data  
sevenyr <- fread("sevenyr_foranalysis_notsc.csv")
sevenyr <- sevenyr %>% select(ID, avg_7yrs, Pollutant)
sixmnth <- fread("6month_avg_exp_foranalysis_ntr.csv")
sixmnth <- sixmnth %>% pivot_longer(!ID, names_to = "Pollutant", values_to = "avg_6mnth")
annual_df3 <- annual_df2 %>% select(ID, Pollutant, days_365_avg)

df <- left_join(sixmnth, annual_df3, by = c("ID", "Pollutant"))
df <- left_join(df, sevenyr, by = c("ID", "Pollutant"))
df2 <- df %>% pivot_longer(!c(ID, Pollutant), names_to = "duration", values_to = "exposure")

df2$duration <- factor(df2$duration, levels = c("avg_6mnth", "days_365_avg", "avg_7yrs"))
str(df2)

p <- df2 %>%
  ggplot(aes(x = duration, y = exposure, fill = Pollutant)) +
  geom_violin(trim = TRUE, draw_quantiles =  c(0.25, 0.5, 0.75), alpha = 0.8) +
  scale_fill_paletteer_d("MetBrewer::Tam") +
  scale_x_discrete(labels = c("6 month", "1 year", "7 years")) +
  # geom_jitter(colour = "black", size = 0.2, alpha = 0.3) +
  theme_minimal() +
  labs(y = "Average exposure over varying time frames", x = "Duration") +
  facet_wrap(~Pollutant, scales = "free_y", nrow = 2) +
  theme(legend.position = "bottom")
ggsave(p, file = "exp_summaries_vplots.pdf", width = 30, height = 15, units = "cm")


#Create a summary table
library(arsenal)
  
annual_df2 <- annual_df2 %>%
  mutate(y_appt = gsub("*....\\|", "", year_range))  
unique(annual_df2$year_range)


control <- tableby.control(,    # Use ANOVA for numeric tests
  total = TRUE,              # Include a total column
  numerics = list(
    mean = FALSE,
    Range = FALSE,
    median = TRUE,           # Show median
    q1q3 = TRUE)) 

tab1 <- summary(tableby(y_appt ~ days_365_avg, data = annual_df2, strata = Pollutant, control = control), text = TRUE)
fwrite(tab1, file = "summary_table.csv")


summary_data <- annual_df2 %>%
  group_by(y_appt, Pollutant) %>%
  summarise(
    count = n(),
    median_value = round(median(days_365_avg, na.rm = TRUE), 2),
    IQR_Value = round(IQR(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(median_value, "[", IQR_Value, "]")
  )
summary_by_pollutant <- annual_df2 %>%
  group_by(Pollutant) %>%
  summarise(
    count = n(),
    median_value = round(median(days_365_avg, na.rm = TRUE), 2),
    IQR_Value = round(IQR(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(median_value, "[", IQR_Value, "]")
  )

summary_wide <- summary_data %>%
  select(-median_value, -IQR_Value, -count) %>%
  pivot_wider(names_from = y_appt, values_from = summary_text, values_fill = "")

fwrite(summary_wide, file = "summary_table2.csv" )


#Also create mean/S.d. summary table
library(arsenal)
  
annual_df2 <- annual_df2 %>%
  mutate(y_appt = gsub("*....\\|", "", year_range))  
unique(annual_df2$year_range)


control <- tableby.control(,    # Use ANOVA for numeric tests
  total = TRUE,              # Include a total column
  numerics = list(
    meansd = TRUE)) 

tab1 <- summary(tableby(y_appt ~ days_365_avg, data = annual_df2, strata = Pollutant, control = control), text = TRUE)
fwrite(tab1, file = "summary_table_meansd.csv")


summary_data <- annual_df2 %>%
  group_by(y_appt, Pollutant) %>%
  summarise(
    count = n(),
    mean_value = round(mean(days_365_avg, na.rm = TRUE), 2),
   sd_value = round(sd(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(mean_value, "[", sd_value, "]")
  )
summary_by_pollutant <- annual_df2 %>%
  group_by(Pollutant) %>%
  summarise(
    count = n(),
    mean_value = round(mean(days_365_avg, na.rm = TRUE), 2),
    sd_Value = round(sd(days_365_avg, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    summary_text = str_c(mean_value, "[", sd_Value, "]")
  )

fwrite(summary_by_pollutant, file = "pollutant_summary_meansd.csv")

