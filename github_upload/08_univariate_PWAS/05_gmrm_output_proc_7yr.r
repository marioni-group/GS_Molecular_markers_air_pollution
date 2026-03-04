#BayesR output processing parallel

library(tidyverse)
library(data.table)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(biomaRt)
library(UpSetR)
library(RColorBrewer)

#Step1_______________________________________________

#Read in step 1 output files & remove test or repeated files
file_list <- list.files(path = "/7yr/output_7yr", recursive = TRUE, pattern = "*.csv", full.names = TRUE)
results_list <- list()

#extract data from files in file-list & create a list of dataframes
for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- fread(file)
  results_list[[filename]] <- data
}
head(results_list)

#Create dataframe of all results
results_df <- bind_rows(results_list, .id = "pollutant")
head(results_df)

#rename output columns as per J.B.
results_df <- results_df %>% 
  dplyr::rename(Iteration = V1,
         Trait_indicator = V2,
         Variance_explained = V3, #?equivalent to SigmaG
         Residual_variance = V4, #equivalent to SigmaE
         Heritability = V5, #This is equal to variance explained/(variance_explained + residual variance) -> %
         No_incl_markers = V6,
         Mixture_group_indicator = V7,
         No_prior_mixture_comp = V8,
         prob_prior_group_1 = V9,
         prob_prior_group_2 = V10,
         prob_prior_group_3 = V11,
         prob_prior_group_4 = V12)
head(results_df)

#Calculate mean heritability on 1250 iterations
heritability <- results_df %>% 
  dplyr::filter(., Iteration >750) %>%
  dplyr::select(Heritability, pollutant) %>%
  dplyr::group_by(pollutant) %>%
  summarise(mean_heritability = mean(Heritability),
          sd_heritability = sd(Heritability))  

#write out heritability results          
fwrite(heritability, file = "/7yr/pollutant_herit_burnin.csv")

#write out full results
fwrite(results_df, file = "/7yr/step1_full_results.csv")


#Step2_______________________________________________________________________________________________

#mlma file =  effect sizes, pip etc
file_list <- list.files(path = "/output_7yr", recursive = TRUE, pattern = "*.mlma", full.names = TRUE)
results_list <- list()

for (file in file_list) {
  filename <- tools::file_path_sans_ext(basename(file))
  data <- fread(file)
  results_list[[filename]] <- data
}
head(results_list)


results_list <- purrr::map(results_list, ~dplyr::rename(., marker_no = V1,
                effect_size = V2,
                pip = V3))

#add in proteins
protein_order <- fread("protein_order.csv")
protein_order <- protein_order %>%
  mutate(marker_no = c(0:132)) #marker should start from 0 
results_list <- purrr::map(results_list, ~ dplyr::left_join(., protein_order, by = "marker_no"))
results_filtered <- purrr::map(results_list, ~dplyr::filter(., pip >0.95))#Select lead CpGs with PIP > 0.95

results_df <- bind_rows(results_list, .id = "pollutant")
head(results_df)
fwrite(results_df, file = "/7yr/full_results_step2.csv")

#Create df for lead CpGs
highpip <- bind_rows(results_filtered, .id = "pollutant") #40 proteins
protein_counts <- highpip %>% dplyr::count(pollutant, sort = TRUE)
pollutant_counts <- highpip %>% dplyr::count(proteins, sort = TRUE)

#write these out:
fwrite(protein_counts, file = "/7yr/Cpg_counts.csv")
fwrite(pollutant_counts, file = "/7yr/pollutant_counts.csv")
fwrite(highpip, file ="/7yr/highpipcpgs.csv")


summary_df <- highpip %>%
  group_by(proteins) %>%
  summarise(associated_traits = paste(pollutant, collapse = ","),
            trait_count = n())
summary_df <- summary_df %>% arrange(desc(trait_count))
summary_df
fwrite(summary_df, file = "/7yr/summary_df.csv")
