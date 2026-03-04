#Summary age_acceleration & pollution work
library(tidyverse)
library(data.table)
library(paletteer)
library(forcats)


#Combine and summarise mixed model results for age-accel - pollution associations

library(tidyverse)
library(data.table)

#_____________________________________________
#Create combined data for 7yr, 6month data and 365day data
#_____________________________________________

#1yr data
poll365 <- fread("/lmekin_models_2ndgen_output.csv")

#6 month data
poll183 <- fread("lmekin_models_output_183poll.csv")
poll183 <- poll183 %>%
  filter(predictor == "PhenoAge_accel" | predictor == "DunedinPACE_accel" | predictor == "DNAmGrimAgev2_accel") %>%
  mutate(data_type = "183_day")

#7yr data
poll7yr <- fread("/lmekin_models_2ndgen_7yrpoll.csv")


#Combine and save out
dfc <- rbind(poll183, poll365, poll7yr)
dfc <- dfc %>% filter(model == "model_4") %>% filter(term == predictor)
fwrite(dfc, file = "model4_diffduratpoll_results.csv")



#First generate plots for 365 data - main text
comdf2 <- fread("lmekin_models_2ndgen_output.csv")

#For supplement
supp <- comdf2 %>% 
  filter(predictor == term) %>% 
  mutate(significance2 = ifelse(p < 0.05/6, "P<0.0083",
                          ifelse(p < 0.05 & p > 0.05/6, "P<0.05", "P>=0.05")))
fwrite(supp, file = "2ndgenclock_results_forsupp_161225.csv")

#Plot final model results
df <- comdf2 %>% 
  filter(predictor == term & model == "model_4") %>%
  mutate(significance2 = ifelse(p < 0.05/6, "P<0.0083",
                          ifelse(p < 0.05 & p > 0.05/6, "P<0.05", "P>=0.05")))


df$predictor <- fct_reorder(df$predictor, df$beta)
df$pollutant <- fct_reorder(df$pollutant, df$beta)

library(scales)
library(paletteer)

p1 <- df %>%
  ggplot(aes(y = pollutant, x = beta, xmin = confint_high, xmax = confint_low, 
             colour = predictor, shape = significance2, position = predictor)) +
  geom_point(position = position_dodge(width = 0.60), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.60), height = .1) +
  scale_color_paletteer_d("MoMAColors::Palermo", labels = c("GrimAge v.2", "DunedinPACE", "PhenoAge")) + 
  scale_shape_manual(values = c(17, 16, 1)) +  # Adjust shape values as needed
  coord_cartesian(xlim = c(-0.06, 0.06)) +
  scale_x_continuous(breaks = seq(-0.06, 0.06, by = 0.02)) +
  labs(title = 'Pollutant exposure and age-acceleration',
       x = 'Effect \nsize', 
       y = 'Pollutant', x = "Beta", colour = "Epigenetic clock \n(age accel)",
       shape = 'Significance') +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', alpha = .5) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white", size = .4),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.title = element_text(size = 14),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 10)),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        strip.placement = "outside") +
        guides(colour = guide_legend(reverse = TRUE))  +
        coord_flip() 

ggsave(p1, file = "model4_plot_100825_flip2_forpresentation.pdf", width = 30, height = 15, units = "cm") 


#_______________________________________________________
#Plots for supplement

#Plot fully adjusted model at different averaging times
df <- fread("model4_diffduratpoll_results.csv") #created above
df <- df %>%
  mutate(significance2 = ifelse(p < 0.05/6, "P<0.0083",
                              ifelse(p < 0.05/2, "P<0.025", 
                          ifelse(p < 0.05 & p > 0.05/2, "P<0.05", "P>=0.05"))))

df$pollutant <- fct_reorder(df$pollutant, df$beta)
df$predictor <- fct_reorder(df$predictor, df$beta)

combplot <- ggplot(data=df, aes(y=pollutant, x=beta, xmin=confint_low, xmax=confint_high, colour = data_type, shape = significance2)) +
 geom_point(position = position_dodge(width = 0.6), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.6), height=.1) +
   scale_color_paletteer_d("NatParksPalettes::Glacier", labels = c("6 months", "1 year", "7 years")) +
   scale_shape_manual(values = c(17, 18, 16, 1)) +
  guides(
    colour = guide_legend("Model")) +
  labs(title='', 
      x='Effect size [95% Confidence Interval]', 
      y = 'Pollutant',
      colour = "Exposure duration",
      shape = "Significance") +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  theme_minimal() + 
  theme(legend.position = "bottom",
        legend.justification = c(0,0),
      #   legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(vjust = 0.2)) +
  facet_wrap(~ predictor, ncol = 3,
             labeller = as_labeller(c(DNAmGrimAgev2_accel = "GrimAge v.2",
                                    DunedinPACE_accel = "DunedinPACE",
                                    PhenoAge_accel = "PhenoAge"))) +
  guides(colour = guide_legend(reverse = TRUE, ncol = 1, title.position = "top"),
      shape = guide_legend(ncol = 1, title.position = "top")) 
ggsave(combplot, file =  "clock_poll_2ndgen_durationcomp_211225.pdf", height = 20, width = 30, units = "cm")


#Plot all models for 365 data
df <- comdf2 %>% 
  filter(predictor == term & !model == "basic") %>%
  dplyr::filter(predictor %in% c("DNAmGrimAgev2_accel", "PhenoAge_accel", "DunedinPACE_accel")) %>%
  mutate(significance2 = ifelse(p < 0.05/6, "P<0.0083",
                              ifelse(p < 0.05/2, "P<0.025", 
                          ifelse(p < 0.05 & p > 0.05/2, "P<0.05", "P>=0.05"))))


df$pollutant <- fct_reorder(df$pollutant, df$beta)

modelplot <- ggplot(data=df, aes(y=pollutant, x=beta, xmin=confint_low, xmax=confint_high, colour = model, shape = significance2)) +
 geom_point(position = position_dodge(width = 0.7), size = 2) + 
  geom_errorbarh(position = position_dodge(width = 0.7), height=.1) +
   scale_color_paletteer_d("PNWColors::Sunset2", labels = c("Model 1", "Model 2", "Model 3", "Model 4")) +
   scale_shape_manual(values = c(17, 18, 16, 1)) +
  guides(
    colour = guide_legend("Model")) +
  labs(title='', x='Effect size [95% Confidence Interval]', y = 'Pollutant',
        colour = "Model",
      shape = "Significance") +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  theme_minimal() + 
  theme(legend.position = "bottom",
         legend.justification = c(0,0),
      #   legend.background = element_rect(fill="white", size=.4),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(vjust = 0.2)) +
  facet_wrap(~ predictor, ncol = 3,
  labeller = as_labeller(c(DNAmGrimAgev2_accel = "GrimAge v.2",
                                    DunedinPACE_accel = "DunedinPACE",
                                    PhenoAge_accel = "PhenoAge"))) +
    guides(colour = guide_legend(ncol = 1, title.position = "top"),
      shape = guide_legend(ncol = 1, title.position = "top")) 
ggsave(modelplot, file = "clock_poll_2ndgen_allmodels_211225.pdf", height = 30, width = 30, units = "cm")
