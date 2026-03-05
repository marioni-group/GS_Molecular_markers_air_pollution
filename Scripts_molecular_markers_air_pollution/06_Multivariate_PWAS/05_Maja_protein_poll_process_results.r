#MAJA_processing

#Maja protein-pollution process results - comparing rintscaled vs scaled pollution

library(tidyverse)
library(data.table)
library(utils)
library(here)


data_path <- ""
betas<- fread(cmd = "unzip -cqp /mean_beta.csv.zip")
betas <- betas[-1,]
prob <- fread("mean_prob.txt")
prob %>% dplyr::filter(V1 > 0.95) %>% dim() 




#add in proteins
proteins <- fread( "proteins_133prepped.csv")
protein_list <- colnames(proteins)
protein_list <- protein_list[-1]
protein_list <- as.data.frame(protein_list)
names(protein_list) <- "protein"
probs <- cbind(protein_list, prob)
probs <- probs %>% dplyr::rename(pip = V1)
highpip <- probs %>% dplyr::filter(pip > 0.95) #46

#add proteins to betas
betas <- cbind(protein_list, betas)

#variance
variance <- fread(cmd = "unzip -cqp var_beta.csv.zip")
head(variance)
variance <- variance[-1,]
variance <- cbind(protein_list, variance)

traits <- fread("annual_avg_exp_foranalysis_scaled.csv") #for order pollutants not IDs
trait_list <- colnames(traits)[-1]
trait_list <- c(trait_list, "logsmok")



#Create results for all input/output
results_full <- left_join(betas, probs, by = c("protein"))
results_full <- results_full %>%
  dplyr::select(protein, pip, everything())

colnames(results_full)[3:11] <- trait_list

results_full2 <- results_full %>%
  pivot_longer(!c(protein, pip),
               names_to = "Trait",
               values_to = "Betas")


colnames(variance)[2:10] <- trait_list
variance <- variance %>%
  pivot_longer(
    !protein,
    names_to = "Trait",
    values_to = "Variance"
  )

results_full2 <- left_join(results_full2, variance, by = c("protein", "Trait"))

results_full2 <- results_full2 %>%
    mutate(stdev = sqrt(Variance))

results_full2 <- results_full2 %>%
    mutate(lower = Betas - stdev,
        upper = Betas + stdev,
        association = ifelse(pip > 0.95 & lower < 0 & upper <0, "yes",
                        ifelse(pip > 0.95 & lower > 0 & upper > 0, "yes", "no")))

fwrite(results_full2, "maja_fullresults_scpoll.csv")

results_full3 <- results_full2 %>%
    mutate(lower2 = Betas - 2*stdev,
        upper2 = Betas + 2*stdev,
        association2 = ifelse(pip > 0.95 & lower2 < 0 & upper2 <0, "yes",
                        ifelse(pip > 0.95 & lower2 > 0 & upper2 > 0, "yes", "no")))
results_full3 %>% filter(association2 == "yes") 
fwrite(results_full3, file = "maja_fullresults_adjthreshold.csv")

summary_df <- results_full3 %>%
  filter(association == "yes") %>%
  group_by(protein) %>%
  summarise(associated_traits = paste(Trait, collapse = ","),
            trait_count = n())

#Add in protein names  
short_annots <- fread("groups_idmapping_01_25_labelled_simplified.csv") #
short_annots <- short_annots %>% rename(protein = original_id)
short_annots <- short_annots %>% dplyr::filter(protein %in% results_plot$protein)
short_annots2 <- short_annots %>% dplyr::select(protein, label)            

summary_df <- left_join(short_annots2, summary_df, by = "protein")


summary_df <- summary_df %>% arrange(desc(trait_count))
fwrite(summary_df, file = "maja_summary_scpoll.csv")



