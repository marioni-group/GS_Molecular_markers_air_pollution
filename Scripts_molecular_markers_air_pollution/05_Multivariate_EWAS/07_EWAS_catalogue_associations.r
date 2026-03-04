#EWAS catalogue summary associations

library(tidyverse)
library(data.table)


#EWAS catalogue search

#1. load ewas catalogue download
#2. load bayesr results
#3. create df of paired CpG and trait from bayesr results 
    #For ewas catalogue:
        #Filter to epigenome-wide significance
        #Filter to whole blood only
        #Have not filtered to N participants as so few protein EWASs
#4. Filter ewas catalogue to CpGs of interest from BayesR
#5. Manually inspect filtered ewas catalogue for traits -> look for different ways of describing proteins
#6. If required:  manually filter ewas catalogue results by trait
#7. And/or create replicated association df using R

#load libraries
library(data.table)
library(R.utils)
library(stringr)

#1 & 2. load bayesr results & ewas catalogue
all_anno <- fread("/maja_highpipresults_1307.csv")
all_anno <- all_anno %>% filter(association == "yes")
catalogue <- fread("ewascatalog-results_050825.txt.gz")
traits <- fread("ewascatalog-studies_050825.txt.gz")

#3 & 4
#Subset to CpGs from bayesR results & to results with epigenome-wide significance
br_cat <- catalogue %>%
  dplyr::filter(CpG %in% all_anno$CpG & P < 3.6e-8)
dim(br_cat)
length(unique(br_cat$CpG)) 
length(unique(br_cat$StudyID)) 

#subset traits file to whole blood and studyIDs in the filtered Cpg results
br_traits <- traits %>%
  dplyr::filter(StudyID %in% br_cat$StudyID & Tissue == "Whole blood")
br_traits <- br_traits %>% filter(N > 1000)
dim(br_traits) 

#merge to common results fulfilling all criteria
both <- merge(br_cat, br_traits, by = "StudyID") 
dim(both) 
length(unique(both$StudyID)) 
length(unique(both$CpG)) 
fwrite(both, file = "ewas_cat_rawres031025.csv")

#summarise traits by CpG
result <- both %>%
  group_by(CpG) %>%
  dplyr::summarise(
    Traits = paste(unique(Trait), collapse = ", "),  # List traits as a comma-separated string
    Num_Traits = n_distinct(Trait))

result <- result %>% arrange(desc(Num_Traits))
fwrite(result, file = "summary_ewascatres031025.csv")


#df of just associations for merging
associations <- all_anno %>%
  dplyr::select(CpG, Trait)

assoc <- associations %>%
  group_by(CpG) %>%
  summarise(
    pollutants = paste(unique(Trait), collapse = ", "),
    number_assoc_pollutants = length(unique(Trait))) %>%
    ungroup() %>%
    distinct()

sumdf <- left_join(assoc, result, by = "CpG")
sumdf <- sumdf %>% arrange(desc(number_assoc_pollutants))
sumdf <- sumdf %>% select(CpG, probe_gene, everything())

#write out both sets of results
fwrite(sumdf, file = "summary_trait_ewascatres031025.csv")

