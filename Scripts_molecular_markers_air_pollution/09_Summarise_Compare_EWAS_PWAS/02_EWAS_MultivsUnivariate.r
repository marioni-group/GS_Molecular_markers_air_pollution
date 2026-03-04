#Compare Multivariate (MAJA) EWAS and Univariate (GMRM) EWAS for 365-day data

library(tidyverse)
library(data.table)
library(qqman)
library(stringr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


results <- fread("/majaewas_highpip_anno_0812_newthresh.csv")


#check cpgs are same for CI thresholds
table(results$association, results$association2)
cpgs1 <- results %>% filter(association == "yes") %>% select(CpG) %>% unique()
cpgs2 <- results %>% filter(association2 == "yes") %>% select(CpG) %>% unique()
identical(cpgs1, cpgs2) #TRUE

#Compare GMRM and MAJA outputs

head(results)
maja <- results
gmrm <- fread("/full_results_step2.csv")

gmrm <- gmrm %>%
  dplyr::rename(Trait = pollutant,
         Betas_gmrm = effect_size,
         pip_gmrm = pip,
         cpg = cpgs)
gmrm_sig <- gmrm %>% 
   filter(pip_gmrm > 0.95) %>%
   mutate(category = ifelse(Trait == "logsmok", "Smoking", "Pollution"))
gmrm_sig <- gmrm_sig %>% mutate(model = "gmrm_1yr")   
fwrite(gmrm_sig, file = "/gmrm_1yr_sigcpgs.csv")

maja <- maja %>%
  dplyr::rename(pip_maja = pip,
        Betas_maja = Betas,
        cpg = CpG)
maja_sig <- maja %>% 
   filter(pip_maja > 0.95 & association2 == "yes") %>%
   mutate(category = ifelse(Trait == "logsmok", "Smoking", "Pollution"))

both <- left_join(gmrm, maja, by = c("Trait", "cpg"))
bothf <- both %>% 
  filter(pip_gmrm > 0.95 | pip_maja > 0.95 & association2 == "yes") #This contains associations for both gmrm and maja
bothf <- bothf %>%
   mutate(maja_association = association,
                gmrm_association = ifelse(pip_gmrm > 0.95, "yes", "no"))  
table(bothf$maja_association)
bothf <- bothf %>% mutate(maja_association = ifelse(is.na(maja_association), "no", maja_association))
table(bothf$maja_association)
table(bothf$gmrm_association)
x <- bothf %>% filter(maja_association == "yes" & gmrm_association == "yes")


#Correlate effect sizes for overlapping associations
cor.test(x$Betas_gmrm, x$Betas_maja) 

#         Pearson's product-moment correlation

# data:  x$Betas_gmrm and x$Betas_maja
# t = 51.534, df = 13, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9924585 0.9992129
# sample estimates:
#       cor 
# 0.9975615 

#save out filtered overlap dataframe
fwrite(bothf, file = "/gmrmmajaewas_comp_highpip_0812_newthresh.csv")


#_____________________________________________________________
#____________________________________________________________________
#Scatter plot
library(paletteer)
clist <- paletteer_d("MetBrewer::Tam")
clist <- c( "#910000FF", clist)
clist

pdf("/scatter_maja_gmrm_cpg_overlap_0912.pdf", family = "Helvetica")
ggplot(bothf, aes(x = Betas_maja, y = Betas_gmrm)) +
  geom_point(aes(color = Trait), size = 3, alpha = 0.8) +
  scale_colour_manual(values = clist) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_minimal() +
  labs(
    x = "Effect Size (MAJA)",
    y = "Effect Size (GMRM)",
    title = "Effect Size Alignment for Overlapping Associations"
  )
dev.off()

#Also make scatter for only overlapping associations:
bothf2 <- bothf %>% filter(maja_association == "yes" & gmrm_association == "yes")
pdf("/scatter_maja_gmrm_cpg_overlap_0912_filt.pdf", family = "Helvetica")
ggplot(bothf2, aes(x = Betas_maja, y = Betas_gmrm)) +
  geom_point(aes(color = Trait), size = 3, alpha = 0.8) +
  scale_colour_manual(values = clist) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_minimal() +
  labs(
    x = "Effect Size (MAJA)",
    y = "Effect Size (GMRM)",
    title = "Effect Size Alignment for Overlapping Associations"
  )
dev.off()



# #____________________________________________________________________________________________________________
# #Add in 6month and 7yr gmrm ewases to results -> compare associations & effect sizes across methods/durations
        
#save out gmrm results for supplement
sixmonth <- fread("/6month/6month_anno_highpipcpgs.csv")
sixdf <- sixmonth %>%
  mutate(exp_dur = "sixmonth") %>%
  dplyr::select(-c(marker_no)) %>%
  mutate(pollutant = str_remove_all(pollutant, "6month"))

oneyr <- fread("/anno_highpipcpgs.csv")
onedf <- oneyr %>%
  mutate(exp_dur = "oneyr") %>%
    dplyr::select(-c(marker_no))

sevenyr <- fread("/7yr/7yr_anno_highpipcpgs.csv")
sevendf <- sevenyr %>%
  mutate(exp_dur = "sevenyr") %>%
  dplyr::select(-c(marker_no)) %>%
  mutate(pollutant = str_remove_all(pollutant, "7yr"))

#Save out for supplementary
gmrm_results <- rbind(onedf, sixdf, sevendf)
dim(gmrm_results)
gmrm_results <- gmrm_results %>% select(exp_dur, pollutant, CpGs, pip, effect_size, chr, UCSC_RefGene_Name, Relation_to_Island)
fwrite(gmrm_results, file = "/summary_results_forsuppl_091225.csv")

#set-up for comparison with MAJA results
sixmonth <- sixmonth %>%
  mutate(exp_dur = "sixmonth") %>%
  dplyr::select(-c(marker_no, pip, strand, chr, pos, Relation_to_Island)) %>%
  mutate(pollutant = str_remove_all(pollutant, "6month"))

oneyr <- oneyr %>%
  mutate(exp_dur = "oneyr") %>%
    dplyr::select(-c(marker_no, pip, strand, chr, pos, Relation_to_Island))


sevenyr <- sevenyr %>%
  mutate(exp_dur = "sevenyr") %>%
  dplyr::select(-c(marker_no, pip, strand, chr, pos, Relation_to_Island)) %>%
  mutate(pollutant = str_remove_all(pollutant, "7yr"))

#Add in MAJA results
maja_1yr <- fread("/majaewas_highpip_anno_0812_newthresh.csv")
maja_1yr <- maja_1yr %>% filter(association2 == "yes")
maja_1yr <- maja_1yr %>%
  mutate(exp_dur = "maja_oneyr")
maja_1yr <- maja_1yr %>% 
    dplyr::select(Trait, Betas, CpG, UCSC_RefGene_Name, exp_dur) %>%
    dplyr::rename(pollutant = Trait,
                  effect_size = Betas,
                  CpGs = CpG)


#Combine all together
df <- rbind(sixmonth, oneyr, sevenyr, maja_1yr)
df2 <- df %>% 
  pivot_wider(
    names_from = exp_dur,
    values_from = effect_size
  )

#______________________________________________________________
#intersect - sum numbers
#________________________________________________________________

head(df2)
df3 <- df2 %>% pivot_longer(cols = c("sixmonth", "oneyr", "sevenyr", "maja_oneyr"), 
                            names_to = "duration",
                            values_to = "effect_size")
df3 <- df3 %>% mutate(association = ifelse(!is.na(effect_size), "yes", "no"))
df4 <- df3 %>% filter(!association == "no") 

res <- df4 %>% group_by(pollutant, CpGs) %>%
  summarise(
    assoc_durations = paste(list(unique(duration)), sep = ","), .groups = "drop"
  )

res %>% filter(assoc_durations == "c(\"oneyr\", \"maja_oneyr\")") %>% dim()


#Restructure df2 for table in supplementary/paper
head(df2)
df2 <- df2 %>% select(pollutant, CpGs, UCSC_RefGene_Name, maja_oneyr, oneyr, sixmonth, sevenyr)
df2 <- df2 %>%
  mutate(n_assocs = rowSums(!is.na(across(where(is.numeric)))))
df2 <- df2 %>% arrange(desc(n_assocs), pollutant, CpGs)

fwrite(df2, file = "/maja_gmrm_overlaptable_newthresh_091225.csv")


#______________________________________________________________
##compare effect sizes - heatmap where fill colour reflects effect size
#________________________________________________________________

#pivot longer for plot
df3 <- df2 %>% pivot_longer(
    cols = c(maja_oneyr, oneyr, sixmonth, sevenyr),
    names_to = "exp_dur",
    values_to = "effect_size"
)

#remove NA values for specific pollutants
df4 <- df3 %>%
  group_by(pollutant, CpGs) %>%
  filter(!is.na(effect_size)) %>%
  ungroup()

library(egg)
p1 <- df4 %>%
    filter(pollutant == "logsmok") %>%
    ggplot(aes(x = exp_dur, y = CpGs, fill = effect_size)) +
   geom_tile(width = 1, height = 0.9) +
   theme_minimal() +
   facet_wrap(~ pollutant, ncol = 1, scales = "free_y", strip.position = "left") +
   theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
   )
p1


p2 <- df4 %>%
    filter(pollutant == "NO3_C") %>%
    ggplot(aes(x = exp_dur, y = CpGs, fill = effect_size)) +
   geom_tile(width = 1, height = 0.9) +
   theme_minimal() +
   facet_wrap(~ pollutant, ncol = 1, scales = "free_y", strip.position = "left") +
   theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0))
p2

# Ensure 'duration' and 'marker' are factors with desired levels
all_durations <- unique(df4$exp_dur)
all_durations <- factor(all_durations, levels = c("sevenyr", "sixmonth", "oneyr", "maja_oneyr"))
df4 <- df4 %>%
  mutate(exp_dur = factor(exp_dur, levels = all_durations))

df5 <- df4 %>%
  group_by(pollutant, CpGs) %>%
  complete(exp_dur = all_durations, fill = list(effect_size = NA)) %>%
  ungroup()


clist <- paletteer_d("rcartocolor::Magenta")
clist2 <- paletteer_d("khroma::BuRd")
clist2 <- clist2[1:8]
clist3 <- paletteer_d("colorBlindness::Blue2Orange12Steps")

fullplot <- ggplot(df5, aes(x = CpGs, y = exp_dur, fill = effect_size)) +
  geom_tile(color = "white", width = 1, height = 1) +  # Fixed tile size
  scale_fill_gradient2(
  low = "#2166ACFF",
  mid = "white",
  high = "#FF6619FF",
  midpoint = 0,
  na.value = "grey90"
) +
  scale_y_discrete(
    labels = c("maja_oneyr" = "MAJA 1yr",
               "oneyr" = "GMRM 1yr",
               "sixmonth" = "GMRM 6mnth",
               "sevenyr" = "GMRM 7yr")
  ) +
  facet_grid(cols = vars(pollutant), 
                scales = "free_x", 
                space = "free_x", 
                switch = "x") +          # Facet by variable with marker axis free
  theme_minimal(base_size = 10) +
  theme(
    # aspect.ratio = 1,                # Make tiles square
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.9, size = 6),
     strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  # coord_fixed(ratio = 1) +
  labs(y = "Exposure Duration", x = "", fill = "Effect Size")
ggsave(fullplot, file = "/univariate_comps2_newthresh.pdf", height = 8, width = 25, units = "cm")


#___________________________________________________________________________
#compare variance explained

maja <- fread("/ewas_variances_0725.csv") 
maja <- maja %>%
  mutate(model = "maja_oneyr") %>%
  dplyr::rename(pollutant = Measure) %>%
  dplyr::select(-c(Variance_of_variance, lower, upper))

sixmonth <- fread("/6month_pollutant_herit_burnin.csv")
sixmonth <- sixmonth %>%
  mutate(model ="gmrm_6mnth",
         pollutant = str_remove_all(pollutant, "6month")) %>%
  dplyr::rename(Variance = mean_heritability) %>%
  mutate(Variance = 100*Variance) %>%
  dplyr::select(-c(sd_heritability))

oneyr <- fread("/pollutant_herit_burnin.csv")
oneyr <- oneyr %>%
  mutate(model ="gmrm_oneyr") %>%
  dplyr::rename(Variance = mean_heritability) %>%
  mutate(Variance = 100*Variance) %>%
  dplyr::select(-c(sd_heritability)) 

sevenyr <- fread("/7yr/7yr_pollutant_variance_expl.csv")
sevenyr <- sevenyr %>%
  mutate(model ="gmrm_7yr",
  pollutant = str_remove_all(pollutant, "7yr")) %>%
  dplyr::rename(Variance = mean_heritability) %>%
    mutate(Variance = 100*Variance) %>%
  dplyr::select(-c(sd_heritability))

df <- rbind(maja, oneyr, sixmonth, sevenyr)
df$model <- factor(df$model, levels = c("maja_oneyr", "gmrm_oneyr", "gmrm_6mnth", "gmrm_7yr"))


df2 <- df %>%
  complete(pollutant, model)


library(paletteer)
ggplot(df2, aes(x = pollutant, y = Variance, fill = model)) +
geom_col(position = "dodge") +
scale_fill_paletteer_d("nationalparkcolors::Acadia", labels = c("MAJA 1yr", "GMRM 1yr", "GMRM 6mnth", "GMRM 7yr")) + 
 guides(
    fill = guide_legend("Exposure duration")) +
labs(x = "Pollutant", y = "Variance explained") +
theme_minimal()



#add in significant CpG Numbers for each model
highpip_results2 <- fread("/maja_highpipresults_1307_newthresh.csv")
cpg_maja <- highpip_results2 %>%
  dplyr::filter(association2 == "yes") %>%
  dplyr::count(Trait, sort = TRUE)

cpg_maja <- cpg_maja %>% 
  rename(pollutant = Trait) %>%
  mutate(model = "maja_oneyr")

gmrm1 <- fread("/Cpg_counts.csv")
gmrm1 <- gmrm1 %>% mutate(
  model = "gmrm_oneyr"
)

gmrm6 <- fread("/6month/6month_Cpg_counts.csv")
gmrm6 <- gmrm6 %>% mutate(model = "gmrm_6mnth")

gmrm7 <- fread("/7yr/7yr_highpipcpgs.csv")
gmrm7 <- gmrm7 %>% select(pollutant) %>% mutate(n = 1) %>% mutate(model = "gmrm_7yr")

countdf <- rbind(cpg_maja, gmrm1, gmrm6, gmrm7)
countdf <- countdf %>%
  complete(pollutant, model)
countdf$pollutant <- str_remove_all(countdf$pollutant, c("6month"))
countdf$pollutant <- str_remove_all(countdf$pollutant, c("7yr"))

#plot together with bar chart

df3 <- left_join(df2, countdf, by = c("pollutant", "model"))
df3$model <- factor(df3$model, levels = c("maja_oneyr", "gmrm_oneyr", "gmrm_6mnth", "gmrm_7yr"))
df3 <- df3 %>% mutate(type = "Count")

vplot <- df3 %>%
ggplot() +
 geom_col(aes(x = pollutant, y = Variance, fill = model), 
 position = position_dodge(width = 0.9)) +
scale_fill_paletteer_d("nationalparkcolors::Acadia", labels = c("MAJA 1yr", "GMRM 1yr", "GMRM 6mnth", "GMRM 7yr")) + 
geom_point(aes(x = pollutant, y = n, colour = type, group = model), 
  position = position_dodge(width = 0.9),
  size = 1, show.legend = TRUE) +
scale_colour_manual(name = "", values = c("Count" = "#E04B28FF"), labels = c("N CpGs")) +
scale_y_continuous(
  name = "Variance explained",
  sec.axis = sec_axis(~ ., name = "Number high pip CpGs")
) +
 guides(
    fill = guide_legend("Exposure duration", override.aes = list(shape = NA)),
    colour = guide_legend("")) +
labs(x = "Pollutant") +
theme_minimal()
ggsave(file = "/variance_comps_newthresh_091225.pdf", width = 20, height = 10, units = "cm")