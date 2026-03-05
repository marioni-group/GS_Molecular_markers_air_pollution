#Compare 
#Maja 1yr results
#GMRM 1yr results
#GMRM 6month results
#GMRM 7yr results

#1. Heat map of exposure/method duration and associations by pollutant
#2. Effect size alignment for overlapping associations (1yr data maja & GMRM)
#3. Exposure duration, variance explained & N associated proteins


#___________________________________________________________________________
#Set up data

#Combine/compare 6month, 1yr, 7yr pollution PWASs (univariate)
library(tidyverse)
library(data.table)
library(stringr)


#load data for CpG associations
sixmonth <- fread("/proteins/6month/highpipcpgs.csv")
sixmonth <- sixmonth %>%
  mutate(exp_dur = "sixmonth") %>%
  dplyr::select(-c(marker_no))

oneyr <- fread("/proteins/1yr/highpipcpgs.csv")
oneyr <- oneyr %>%
  mutate(exp_dur = "oneyr") %>%
    dplyr::select(-c(marker_no))
#explore for paper
head(oneyr)    
dim(oneyr) #189
oneyr %>% filter(pollutant == "logsmok") %>% dim()
oneyr %>% count(pollutant, sort = TRUE)

sevenyr <- fread("/proteins/7yr/highpipcpgs.csv")
sevenyr <- sevenyr %>%
  mutate(exp_dur = "sevenyr") %>%
  dplyr::select(-c(marker_no))

#Combine to save out for supplementary & summarise
df <- rbind(sixmonth, oneyr, sevenyr)
df <- df %>% dplyr::select(exp_dur, pollutant, proteins, pip, effect_size)
df %>% group_by(exp_dur) %>% count(pollutant, sort = TRUE)

dfsum <- df %>%
  group_by(exp_dur, pollutant) %>%
  summarise(
    pollutants = paste(unique(proteins), collapse = ", "),
    number_assoc_proteins = length(unique(proteins))) %>%
    ungroup() %>%
    distinct()
fwrite(dfsum, file = "/proteins/gmrm_results_summary.csv")

#add annotations
short_annots <- fread("groups_idmapping_01_25_labelled_simplified.csv") #
short_annots <- short_annots %>% rename(protein = original_id)
short_annots <- short_annots %>% dplyr::filter(protein %in% df$proteins)
df <- df %>% dplyr::rename(protein = proteins)
df2 <- left_join(df, short_annots)
df2 <- df2 %>% select(-c(Entry, From, Reviewed, "Entry Name", protein_names, Organism, Length))
fwrite(df2, file = "/proteins/gmrm_protresults_allexp.csv")



maja_1yr <- fread("/maja_fullresults_scpoll_burninadj_0812_newthresh.csv")
maja_1yr <- maja_1yr %>% filter(association2 == "yes")
maja_1yr <- maja_1yr %>%
  mutate(exp_dur = "maja_oneyr")
maja_1yr <- maja_1yr %>% 
    dplyr::select(Trait, Betas, pip, protein, exp_dur) %>%
    dplyr::rename(pollutant = Trait,
                  effect_size = Betas,
                  proteins = protein)

#Combine
df <- rbind(sixmonth, oneyr, sevenyr, maja_1yr)
df2 <- df %>% 
  pivot_wider(
    names_from = exp_dur,
    values_from = effect_size
  )


#Restructure df2 for table in supplementary/paper
head(df2)
df2 <- df2 %>% select(pollutant, proteins, maja_oneyr, oneyr, sixmonth, sevenyr)
df2 <- df2 %>%
  mutate(n_assocs = rowSums(!is.na(across(where(is.numeric)))))
df2 <- df2 %>% arrange(desc(n_assocs), pollutant, proteins)

#Annotate with protein gene names & labels
short_annots <- fread("groups_idmapping_01_25_labelled_simplified.csv") #
short_annots <- short_annots %>% rename(protein = original_id)
short_annots <- short_annots %>% dplyr::filter(protein %in% df2$proteins)
short_annots2 <- short_annots %>% dplyr::select(protein, label, gene_names_simpl) 
short_annots2 <- short_annots2 %>% rename(proteins = protein)

df3 <- left_join(df2, short_annots2, by = "proteins")
df3 <- df3 %>% select(pollutant, proteins, gene_names_simpl, label, everything())

fwrite(df3, file = "maja_gmrm_overlaptable_newthresh_101225.csv")


#_________________________________________________________________________
#1. Heat map of exposure/method duration and associations by pollutant

#pivot longer for plot
df4 <- df3 %>% pivot_longer(
    cols = c(maja_oneyr, oneyr, sixmonth, sevenyr),
    names_to = "exp_dur",
    values_to = "effect_size"
)

#remove NA values for specific pollutants
df4 <- df4 %>%
  group_by(pollutant, proteins, gene_names_simpl) %>%
  select(-c(n_assocs, label)) %>%
  filter(!is.na(effect_size)) %>%
  ungroup()

  # Ensure 'duration' and 'marker' are factors with desired levels
all_durations <- unique(df4$exp_dur)
all_durations <- factor(all_durations, levels = c("sevenyr", "sixmonth", "oneyr", "maja_oneyr"))
df4 <- df4 %>%
  mutate(exp_dur = factor(exp_dur, levels = all_durations))

df5 <- df4 %>%
  group_by(pollutant, proteins, gene_names_simpl) %>%
  complete(exp_dur = all_durations, fill = list(effect_size = NA)) %>%
  ungroup()

#check 
check <- df5 %>%
  group_by(pollutant, proteins, gene_names_simpl) %>%
  filter(!is.na(effect_size))


fullplot <- ggplot(df5, aes(x = gene_names_simpl, y = exp_dur, fill = effect_size)) +
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
  theme_minimal(base_size = 8) +
  theme(
    # aspect.ratio = 1,                # Make tiles square
    axis.text.x = element_text(angle = 90, hjust = 0.9, vjust = 0.9, size = 2),
     strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.key.size = unit(0.3, "cm")
  ) +
  labs(y = "Exposure Duration", x = "", fill = "Effect Size")  
ggsave("/proteins/maja_gmrm_effectsizes3_newthresh_101225.pdf", width = 30, height = 8, units = "cm")


#Make a vertical plot

#adjust factor level
durations <- factor(all_durations, levels = c("maja_oneyr", "oneyr", "sixmonth", "sevenyr"))

df6 <- df5 %>% mutate(exp_dur = factor(exp_dur, levels = durations))
df6 <- df6 %>%
  group_by(pollutant, proteins, gene_names_simpl) %>%
  complete(exp_dur = all_durations, fill = list(effect_size = NA)) %>%
  ungroup()

ggplot(df5, aes(x = exp_dur, y = gene_names_simpl, fill = effect_size)) +
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
  facet_grid(rows = vars(pollutant), 
                scales = "free_y", 
                space = "free_y", 
                switch = "y") +          # Facet by variable with marker axis free
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 2),
     strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 6),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  guides(fill = guide_colorbar(barheight = unit(25, "mm"),
                               barwidth  = unit(3,  "mm"))) +
theme(legend.title = element_text(size = 7),
       legend.text  = element_text(size = 6),
       legend.box.spacing = unit(1, "mm")) +
  labs(y = "Protein genes", x = "Exposure duration", fill = "Effect Size")  
ggsave("proteins/maja_gmrm_effectsizes4_newthresh.pdf", width = 10, height = 40, units = "cm")


#___________________________________________________________________________
#compare variance explained


maja <- fread("maja_prot_variance_explained.csv") 
maja <- maja %>%
  mutate(model = "maja_oneyr") %>%
  dplyr::rename(pollutant = Measure) %>%
  dplyr::select(-c(Variance_of_variance, lower, upper))

sixmonth <- fread("/6month/pollutant_herit_burnin.csv")
sixmonth <- sixmonth %>%
  mutate(model ="gmrm_6mnth",
         pollutant = str_remove_all(pollutant, "6month")) %>%
  dplyr::rename(Variance = mean_heritability) %>%
  mutate(Variance = 100*Variance) %>%
  dplyr::select(-c(sd_heritability))

oneyr <- fread("/1yr/pollutant_herit_burnin.csv")
oneyr <- oneyr %>%
  mutate(model ="gmrm_oneyr") %>%
  dplyr::rename(Variance = mean_heritability) %>%
  mutate(Variance = 100*Variance) %>%
  dplyr::select(-c(sd_heritability)) 

sevenyr <- fread("/7yr/pollutant_herit_burnin.csv")
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


#fill empties for logsmok
#change colour scheme
#include var of var or credible intervals?
#add N of CPGS in Ewas
library(paletteer)
ggplot(df2, aes(x = pollutant, y = Variance, fill = model)) +
geom_col(position = "dodge") +
scale_fill_paletteer_d("nationalparkcolors::Acadia", labels = c("MAJA 1yr", "GMRM 1yr", "GMRM 6mnth", "GMRM 7yr")) + 
 guides(
    fill = guide_legend("Exposure duration")) +
labs(x = "Pollutant", y = "Variance explained") +
theme_minimal()


#add in significant CpG Numbers for each model
prot_maja <- fread("maja_fullresults_scpoll_burninadj_0812_newthresh.csv")

prot_maja2 <- prot_maja %>%   
  filter(association2 == "yes") %>% 
  dplyr::count(Trait, sort = TRUE) %>%
  rename(pollutant = Trait) %>%
  mutate(model = "maja_oneyr")

gmrm1 <- fread("/1yr/pollutant_counts.csv")
gmrm1 <- gmrm1 %>% mutate(
  model = "gmrm_oneyr"
)

gmrm6 <- fread("/6month/highpipcpgs.csv")
gmrm6 <- gmrm6 %>% 
  dplyr::count(pollutant, sort = TRUE) %>%
  mutate(model = "gmrm_6mnth")

gmrm7 <- fread("/7yr/Cpg_counts.csv")
gmrm7 <- gmrm7 %>% mutate(model = "gmrm_7yr")

countdf <- rbind(prot_maja2, gmrm1, gmrm6, gmrm7)
countdf <- countdf %>%
  complete(pollutant, model)


#plot together with bar chart

df3 <- left_join(df2, countdf, by = c("pollutant", "model"))
df3$model <- factor(df3$model, levels = c("maja_oneyr", "gmrm_oneyr", "gmrm_6mnth", "gmrm_7yr"))
df3 <- df3 %>% mutate(type = "Count")

vplot <- df3 %>%
ggplot() +
 geom_col(aes(x = pollutant, y = Variance, fill = model), 
 position = position_dodge(width = 0.9)) +
scale_fill_paletteer_d("nationalparkcolors::Acadia", labels = c("MAJA 1yr", "GMRM 1yr", "GMRM 6mnth", "GMRM 7yr")) + 
geom_point(aes(x = pollutant, y = n, colour = model, group = model), 
  position = position_dodge(width = 0.9),
  size = 1, show.legend = TRUE) +
scale_colour_paletteer_d("nationalparkcolors::Acadia") +
scale_y_continuous(
  name = "Variance explained",
  sec.axis = sec_axis(~ ., name = "Number high pip CpGs")
) +
 guides(
    fill = guide_legend("Exposure duration", override.aes = list(shape = NA)),
    colour = guide_legend("")) +
labs(x = "Pollutant") +
theme_minimal()
ggsave(vplot, file = "/proteins/prot_variance_comps_newthresh_101225.pdf", width = 20, height = 10, units = "cm")



#____________________________________________________________________
#Scatter plot of 1 year MAJA and 1 year gmrm result effect sizes

#Compare GMRM and MAJA outputs: set up data
head(results)
maja <- fread("/maja_fullresults_scpoll_burninadj_0812_newthresh.csv")
gmrm <- fread("/1yr/full_results_step2.csv")

gmrm <- gmrm %>%
  dplyr::rename(Trait = pollutant,
         Betas_gmrm = effect_size,
         pip_gmrm = pip,
         protein = proteins)
gmrm_sig <- gmrm %>% 
   filter(pip_gmrm > 0.95) %>%
   mutate(category = ifelse(Trait == "logsmok", "Smoking", "Pollution"))
gmrm_sig <- gmrm_sig %>% mutate(model = "gmrm_1yr")   

maja <- maja %>%
  filter(association2 == "yes") %>%
  select(-c(Variance, stdev, lower, upper)) %>%
  dplyr::rename(pip_maja = pip,
        Betas_maja = Betas)
        
maja_sig <- maja %>% 
   filter(pip_maja > 0.95 & association == "yes") %>%
   mutate(category = ifelse(Trait == "logsmok", "Smoking", "Pollution"))

both <- dplyr::full_join(gmrm, maja, by = c("Trait", "protein"))
bothf <- both %>% 
  filter(pip_gmrm > 0.95 | pip_maja > 0.95 & association == "yes") #This contains associations for both gmrm and maja
bothf <- bothf %>%
   mutate(maja_association = association,
                gmrm_association = ifelse(pip_gmrm > 0.95, "yes", "no"),
          )  %>%
    filter(!gmrm_association == "no") 

#run this step if want to include associations unique to each model
bothf2 <- bothf %>%
  mutate(
    Betas_maja = ifelse(is.na(Betas_maja), 0, Betas_maja),
    Betas_gmrm = ifelse(is.na(Betas_gmrm), 0, Betas_gmrm),
    maja_association = ifelse(is.na(maja_association), "no", maja_association)
  ) 

table(bothf2$maja_association)
#  no yes 
#  62 127

table(unique(maja_sig$protein) %in% unique(gmrm_sig$protein)) 

# FALSE  TRUE 
#     2    43 

table(unique(gmrm_sig$protein) %in% unique(maja_sig$protein))
# FALSE  TRUE 
#    36    43 


table(bothf$gmrm_association)
#  no yes 
# 147 189 



library(paletteer)
clist <- paletteer_d("MetBrewer::Tam")
clist <- c( "#910000FF", clist)
clist


ggplot(bothf, aes(x = Betas_maja, y = Betas_gmrm)) +
  geom_point(aes(color = Trait), size = 2, alpha = 0.7) +
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

ggsave("effect_size_comp_overlaponly_updated_101225.pdf", width = 20, height = 20, units = "cm", family = "Helvetica")  

cor(bothf$Betas_maja, bothf$Betas_gmrm, use = "complete.obs")
# [1] 0.9895601


ggplot(bothf2, aes(x = Betas_maja, y = Betas_gmrm)) +
  geom_point(aes(color = Trait), size = 2, alpha = 0.7) +
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

ggsave("/effect_size_comp_updated_101225.pdf", family = "helvetica") 

