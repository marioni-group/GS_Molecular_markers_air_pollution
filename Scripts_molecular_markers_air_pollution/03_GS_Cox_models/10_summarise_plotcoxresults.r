#summarise & plot cvd & HTN cox model results
library(data.table)
library(tidyverse)

#Load results
cvd_df <- fread("coxme_cvdagebmipollcat_0707.csv")
cvd_df <- cvd_df %>% mutate(outcome = "cvd")
htn_df <- fread("htn_coxme_agebmipollcat_0707.csv")
htn_df <- htn_df %>% mutate(outcome = "htn")
mi_df  <- fread("coxme_agebmipollcat_mi_0707.csv")
mi_df <- mi_df %>% mutate(outcome = "mi")
ihd_df <- fread("coxme_agebmipollcat_ihd_0707.csv")
ihd_df <- ihd_df %>% mutate(outcome = "ihd")
stroke_df <- fread("coxme_agebmipollcat_stroke_0707.csv")
stroke_df <- stroke_df %>% mutate(outcome = "stroke")
dementia_df <- fread("/Dementia/cox_results_pollcat_full_0807.csv")
dementia_df <- dementia_df %>% mutate(outcome = "dementia")


#Filter to pollutant terms only
df <- rbind(cvd_df, htn_df, mi_df, ihd_df, stroke_df, dementia_df)
df2 <- df %>% filter(str_detect(term, pollutant))
df2 <- df2 %>% 
  mutate(significance = ifelse(p.value < 0.05 & p.value > 0.05/12, "nominal", 
                            ifelse(p.value < 0.05/12, "bonf_string", "none")))

df2 <- df2 %>% mutate(pollutant = str_sub(pollutant, end = -2),
                      term = str_sub(term, end = -4))
df2 <- df2 %>% select(outcome, everything())
fwrite(df2, file = "cox_all_results_forsupp.csv")


sig <- df2 %>% 
  group_by(outcome) %>% 
  dplyr::filter(significance == "bonf_string" & model == "model_3") %>%
  summarise(
    outcome = outcome,
    pollutant = str_sub(pollutant, end = -2),
    HR = estimate,
    CI = paste0(round(conf.low, 2), "-", round(conf.high, 2)),
    P_value = round(p.value, 3))
fwrite(sig, file = "summary_coxmodels_sig.csv")


#filter to pollutant term only
df <- df %>%
  filter(term == paste0(pollutant, "Q2")) %>%
  mutate(significance = ifelse(p.value < 0.05/12, "P<0.004",
                          ifelse(p.value < 0.05 & p.value > 0.05/12, "P<0.05", "P>=0.05")))

df$pollutant <- str_sub(df$pollutant, end = -2)


#Plot with all models
p1 <- df %>%
ggplot(aes(y=pollutant, x=estimate, xmin=conf.low, xmax=conf.high, colour = model, shape = significance)) +
  geom_point(position = position_dodge(width = 0.75), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.1) +
  scale_color_paletteer_d("PNWColors::Sunset2") +
  scale_shape_manual(values = c(17, 18, 16, 1)) +
    coord_cartesian(xlim = c(0.6, 2.6)) +
  scale_x_continuous(breaks = seq(0.6, 2.6, by = 0.4)) +
  labs(title='',
     x='Hazard Ratio [95% Confidence Interval]', 
     y = 'Outcome', colour = "model") +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_bw() + 
  theme(legend.position = "right",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title = element_text(size = 10),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(vjust = -0.5),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        strip.placement = "outside") +
  facet_wrap(~outcome,
  labeller = labeller(
    outcome = c("htn" = "Hypertension",
                "cvd" = "CVD",
                "stroke" = "Ischaemic Stroke",
                "mi" = "Myocardial Infarction",
                "ihd" = "Coronary Heart Disease",
                "dementia" = "Dementia")
  ))
ggsave(p1, file ="allmodels_summaryplot.pdf", width = 25, height = 30, units = "cm")

#Plot with final model only
df$outcome <- factor(df$outcome, levels = c("dementia", "stroke", "htn", "ihd", "mi", "cvd"))
df <- df %>% rename(pollutant = poll)
df$pollutant  <- fct_reorder(df$pollutant, df$estimate)

#Pick some colours
library(paletteer)

scale_colour_paletteer_d("MetBrewer::Benedictus")
scale_color_paletteer_d("MetBrewer::Benedictus")
scale_fill_paletteer_d("MetBrewer::Benedictus")
colour_list <- paletteer_d("MetBrewer::Benedictus")

col <- paletteer_d("NatParksPalettes::Glacier")
col <- col[1:4]

#add column with case numbers
df <- df %>% 
  mutate(case_lable = ifelse(outcome == "cvd", "Cases:905 \nControls:12118",
                        ifelse(outcome == "htn", "Cases:1,689 \nControls:19064",
                          ifelse(outcome == "mi", "Cases:145 \nControls:12878",
                            ifelse(outcome == "ihd", "Cases:530 \nControls:12493",
                              ifelse(outcome == "stroke", "Cases:159 \nControls:12878", "Cases:323 \nControls:10416"))))))

# Decide where to place the label inside each facet:
# - x position just inside the left x-limit
# - y position at the top (or bottom) pollutant, nudged a bit to avoid overlap
x_left <- 2.2  
poll_levels <- levels(fct_inorder(df$pollutant))
y_bottom <- poll_levels[1] # bottom row 

# One row per facet (outcome) with a single case_lable
labels_df <- df %>%
  group_by(outcome) %>%
  summarise(case_lable = first(na.omit(case_lable)), .groups = "drop") %>%
  mutate(
    estimate  = x_left,
    pollutant = y_bottom,
  )



p2 <- df %>% filter(model == "model_3") %>%
  ggplot(aes(y=pollutant, x=estimate, xmin=conf.low, xmax=conf.high, shape = significance, colour = significance)) +
  geom_point(position = position_dodge(width = 0.75), size = 2.75) + 
  geom_errorbarh(position = position_dodge(width = 0.75), height=.2, linewidth = 1.05) +
  scale_color_manual(name = "Significance", values = col) +
  scale_shape_manual(name = "Significance", values = c(17, 16, 1)) +
    coord_cartesian(xlim = c(0.6, 2.6)) +
  scale_x_continuous(breaks = seq(0.6, 2.6, by = 0.4)) +
    labs(title='',
     x='Hazard Ratio [95% Confidence Interval]', 
     y = 'Outcome', colour = "Significance") +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_bw() + 
  theme(legend.position = "right",
        legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title = element_text(size = 14),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        strip.placement = "outside") +
  facet_wrap(~outcome,
  labeller = labeller(
    outcome = c("htn" = "Hypertension",
                "cvd" = "Cardiovascular Disease",
                "stroke" = "Ischaemic Stroke",
                "mi" = "Myocardial Infarction",
                "ihd" = "Ischaemic Heart Disease",
                "dementia" = "Dementia")
  )) +
    geom_text(
    data = labels_df,
    aes(x = estimate, y = pollutant, label = case_lable),
    inherit.aes = FALSE,
    hjust = 0.4,              # left-align at x_left
    vjust = 9.8,           # nudge a little above the chosen row
    size = 3
  ) +
  theme(strip.text.x = element_text(size = 16))

ggsave(p2, file = "model3_summaryplot2_160126_fortalk.pdf", width = 30, height = 20, units = "cm")




#Create summary of proportional hazards

stroke <- fread("cox_assump_agebmipollcat__stroke_0707.csv")
stroke <- stroke %>%
  mutate(outcome = "stroke",
         term = str_remove(term, "\\.+\\d+$"),
         pollutant = str_sub(pollutant, end = -2)) %>%
  select(outcome, model, pollutant, term, chisq, df, p) %>%
  mutate(term = ifelse(str_detect(term, pollutant), pollutant, term))


ihd <- fread("cox_assump_agebmipollcat_ihd_0707.csv")
ihd <- ihd %>%
  mutate(outcome = "ihd",
         term = str_remove(term, "\\.+\\d+$"),
         pollutant = str_sub(pollutant, end = -2)) %>%
  select(outcome, model, pollutant, term, chisq, df, p) %>%
  mutate(term = ifelse(str_detect(term, pollutant), pollutant, term))

mi <- fread("cox_assump_agebmipollcat_mi_0707.csv")
mi <- mi %>%
  mutate(outcome = "mi",
         term = str_remove(term, "\\.+\\d+$"),
         pollutant = str_sub(pollutant, end = -2)) %>%
  select(outcome, model, pollutant, term, chisq, df, p) %>%
  mutate(term = ifelse(str_detect(term, pollutant), pollutant, term))


htn <- fread("htn_cox_assump_bmipollcat_0807.csv") #non-age filtered data
htn <- htn %>%
  mutate(outcome = "htn",
         term = str_remove(term, "\\.+\\d+$"),
         pollutant = str_sub(pollutant, end = -2)) %>%
  select(outcome, model, pollutant, term, chisq, df, p) %>%
  mutate(term = ifelse(str_detect(term, pollutant), pollutant, term))


dementia <- fread("/Dementia/cox_assump_pollcat_full_0807.csv")
dementia <- dementia %>%
  mutate(outcome = "dementia",
         term = str_remove(term, "\\.+\\d+$"),
         pollutant = str_sub(pollutant, end = -2)) %>%
  select(outcome, model, pollutant, term, chisq, df, p) %>%
  mutate(term = ifelse(str_detect(term, pollutant), pollutant, term))


cvd <- fread("coxme_assump_cvdagebmipollcat_0707.csv")
cvd <- cvd %>%
  mutate(outcome = "cvd",
         term = str_remove(term, "\\.+\\d+$"),
         pollutant = str_sub(pollutant, end = -2)) %>%
  select(outcome, model, pollutant, term, chisq, df, p) %>%
  mutate(term = ifelse(str_detect(term, pollutant), pollutant, term))


#combine
comb <- rbind(cvd, mi, ihd, stroke, htn, dementia)

#filter to term/global results only
comb2 <- comb %>%
  filter(term == pollutant | term == "GLOBAL")  

comb2 %>% filter(p < 0.05)

fwrite(comb2, file = "coxmodelassumpall.csv")
