#Process EpiScore training step

library(tidyverse)
library(data.table)
library(utils)
library(here)


#Load betas
betas<- fread(cmd = "unzip -cqp /mean_beta.csv.zip")
betas <- betas[-1,]
prob <- fread("mean_prob.txt")
prob %>% dplyr::filter(V1 > 0.95) %>% dim() #19

#add in cpgs
cpgs <- fread("cpg_order_mvalinput.csv")
colnames(cpgs) <- "CpG"
probs <- cbind(cpgs, prob)
probs <- probs %>% dplyr::rename(pip = V1)
highpip <- probs %>% dplyr::filter(pip > 0.95)

#add cpgs to betas
betas <- cbind(cpgs, betas)
head(betas)

#Extract pollutant names
poll <- fread("poll_scaled_15516.csv")
trait_list <- colnames(poll)[2:10]

colnames(betas)[2:10] <- trait_list

#Save out
fwrite(betas, file = "cpg_weights.csv" )

#save out only non-zero cpgs for supplement
betas_filt <- betas %>%
  dplyr::filter(!if_all(c(PM25, NO2, NO, O3, PM10, SO2, NO3_C, NO3_F, logsmok), ~ . == 0))
fwrite(betas_filt, file = "poll_epi_nonzerocpgs.csv")


#summarise number of CpGs with non-0 weights per pollutant
betas2 <- betas %>%
  pivot_longer(
    cols = !CpG,
    names_to = "pollutant",
    values_to = "weights"
  )

epi_cpgs <- betas2 %>% 
  group_by(pollutant) %>%
  summarise(
    weighted_cpgs = n(),
   positive_weights = sum(weights > 0),
    negative_weights = sum(weights < 0)) %>%
  ungroup()
fwrite(epi_cpgs, file = "summary_episcore_cpgs.csv")  


