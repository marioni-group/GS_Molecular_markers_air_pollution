#BayesR output processing parallel

library(tidyverse)
library(data.table)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(biomaRt)
library(UpSetR)
library(RColorBrewer)

#_______________
#1-year data#
#_______________

#Step1_______________________________________________

#Read in step 1 output files & remove test or repeated files
file_list <- list.files(path = "", recursive = TRUE, pattern = "*.csv", full.names = TRUE)
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
         Heritability = V5,
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
fwrite(heritability, file = "/pollutant_herit_burnin.csv")
her <- fread("/pollutant_herit_burnin.csv")
#write out full results
fwrite(results_df, file = "/step1_full_results.csv")


#Step2_______________________________________________________________________________________________

#mlma file =  effect sizes, pip etc
file_list <- list.files(path = "", recursive = TRUE, pattern = "*.mlma", full.names = TRUE)
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

#add in CpGs
CpG_order <- fread("cpg_order_methinput.csv")
CpG_order <- CpG_order %>%
  mutate(marker_no = c(0:752721)) #marker should start from 0 
results_list <- purrr::map(results_list, ~ dplyr::left_join(., CpG_order, by = "marker_no"))
results_filtered <- purrr::map(results_list, ~dplyr::filter(., pip >0.95))#Select lead CpGs with PIP > 0.95

results_df <- bind_rows(results_list, .id = "pollutant")
head(results_df)
fwrite(results_df, file = "/full_results_step2.csv")

#Create df for lead CpGs
highpip <- bind_rows(results_filtered, .id = "pollutant")
Cpg_counts <- highpip %>% dplyr::count(pollutant, sort = TRUE)
pollutant_counts <- highpip %>% dplyr::count(cpgs, sort = TRUE)

#write these out:
fwrite(Cpg_counts, file = "Cpg_counts.csv")
fwrite(protein_counts, file = "pollutant_counts.csv")
fwrite(highpip, file ="highpipcpgs.csv")

#annotate
highpip <- fread("highpipcpgs.csv")
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- data.frame(epic_anno)
setDT(anno, keep.rownames = "CpG")
highpip <- highpip %>% dplyr::rename(CpGs = cpgs)
anno <- anno %>% filter(., CpG %in% highpip$CpGs)
anno <- anno %>% dplyr::select(Name, chr, pos, strand, UCSC_RefGene_Name, Relation_to_Island)

anno$chr <- gsub("chr", "", anno$chr) # remove chr from chromosome names
anno$chr <- gsub("X", "23", anno$chr) # replace chr X and Y with 23 and 24
anno$chr <- gsub("Y", "24", anno$chr)
anno$chr <- as.numeric(anno$chr) # set chr as numeric 

anno$UCSC_RefGene_Name	 <- as.factor(anno$UCSC_RefGene_Name) # set genes and strands as factors 
anno$strand <- as.factor(anno$strand)
anno[which(anno$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA #set missing gene names to NA 

anno$CpGs <- anno$Name
anno <- anno %>% dplyr::select(-Name)

cpg_anno <- left_join(highpip, anno, by = "CpGs")
fwrite(cpg_anno, file = "anno_highpipcpgs.csv")


#_______________
#7-year data#
#_______________

#Step1_______________________________________________

#Read in step 1 output files & remove test or repeated files
file_list <- list.files(path = "/output_7yr", recursive = TRUE, pattern = "*.csv", full.names = TRUE)
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
         Heritability = V5,
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

# [1] "NO26month"   "NO3_C6month" "NO3_F6month" "NO6month"    "O36month"   
# [6] "PM106month"  "PM256month"  "SO26month" 



#write out heritability results          
fwrite(heritability, file = "/7yr_pollutant_variance_expl.csv")
her <- fread("/pollutant_herit_burnin.csv")
#write out full results
fwrite(results_df, file = "7yr_step1_full_results.csv")


#Step2_______________________________________________________________________________________________

#mlma file =  effect sizes, pip etc
file_list <- list.files(path = "output_7yr", recursive = TRUE, pattern = "*.mlma", full.names = TRUE)
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

#add in CpGs
CpG_order <- fread("cpg_order_methinput.csv")
CpG_order <- CpG_order %>%
  mutate(marker_no = c(0:752721)) #marker should start from 0 
results_list <- purrr::map(results_list, ~ dplyr::left_join(., CpG_order, by = "marker_no"))
results_filtered <- purrr::map(results_list, ~dplyr::filter(., pip >0.95))#Select lead CpGs with PIP > 0.95

results_df <- bind_rows(results_list, .id = "pollutant")
head(results_df)
fwrite(results_df, file = "/7yr_full_results_step2.csv")

#Create df for lead CpGs
highpip <- bind_rows(results_filtered, .id = "pollutant")


#write these out:
fwrite(highpip, file ="/7yr_highpipcpgs.csv")

#annotate
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- data.frame(epic_anno)
setDT(anno, keep.rownames = "CpG")
highpip <- highpip %>% dplyr::rename(CpGs = cpgs)
anno <- anno %>% filter(., CpG %in% highpip$CpGs)
anno <- anno %>% dplyr::select(Name, chr, pos, strand, UCSC_RefGene_Name, Relation_to_Island)

anno$chr <- gsub("chr", "", anno$chr) # remove chr from chromosome names
anno$chr <- gsub("X", "23", anno$chr) # replace chr X and Y with 23 and 24
anno$chr <- gsub("Y", "24", anno$chr)
anno$chr <- as.numeric(anno$chr) # set chr as numeric 

anno$UCSC_RefGene_Name	 <- as.factor(anno$UCSC_RefGene_Name) # set genes and strands as factors 
anno$strand <- as.factor(anno$strand)
anno[which(anno$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA #set missing gene names to NA 

anno$CpGs <- anno$Name
anno <- anno %>% dplyr::select(-Name)

cpg_anno <- left_join(highpip, anno, by = "CpGs")
fwrite(cpg_anno, file = "/7yr/7yr_anno_highpipcpgs.csv")


#_______________
#6-month data#
#_______________

#Step1_______________________________________________

#Read in step 1 output files & remove test or repeated files
file_list <- list.files(path = "/output_6month", recursive = TRUE, pattern = "*.csv", full.names = TRUE)
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
         Heritability = V5,
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

# [1] "NO26month"   "NO3_C6month" "NO3_F6month" "NO6month"    "O36month"   
# [6] "PM106month"  "PM256month"  "SO26month" 



#write out heritability results          
fwrite(heritability, file = "/6month/6month_pollutant_herit_burnin.csv")

#write out full results
fwrite(results_df, file = "/6month/6month_step1_full_results.csv")


#Step2_______________________________________________________________________________________________

#mlma file =  effect sizes, pip etc
file_list <- list.files(path = "output_6month", recursive = TRUE, pattern = "*.mlma", full.names = TRUE)
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

#add in CpGs
CpG_order <- fread("cpg_order_methinput.csv")
CpG_order <- CpG_order %>%
  mutate(marker_no = c(0:752721)) #marker should start from 0 
results_list <- purrr::map(results_list, ~ dplyr::left_join(., CpG_order, by = "marker_no"))
results_filtered <- purrr::map(results_list, ~dplyr::filter(., pip >0.95))#Select lead CpGs with PIP > 0.95

results_df <- bind_rows(results_list, .id = "pollutant")
head(results_df)
fwrite(results_df, file = "6month_full_results_step2.csv")

#Create df for lead CpGs
highpip <- bind_rows(results_filtered, .id = "pollutant")
Cpg_counts <- highpip %>% dplyr::count(pollutant, sort = TRUE)
pollutant_counts <- highpip %>% dplyr::count(cpgs, sort = TRUE)

#write these out:
fwrite(Cpg_counts, file = "/6month/6month_Cpg_counts.csv")
fwrite(pollutant_counts, file = "/6month/6month_pollutant_counts.csv")
fwrite(highpip, file ="/6month/6month_highpipcpgs.csv")

#annotate
highpip <- fread("/6month/6month_highpipcpgs.csv")
epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- data.frame(epic_anno)
setDT(anno, keep.rownames = "CpG")
highpip <- highpip %>% dplyr::rename(CpGs = cpgs)
anno <- anno %>% filter(., CpG %in% highpip$CpGs)
anno <- anno %>% dplyr::select(Name, chr, pos, strand, UCSC_RefGene_Name, Relation_to_Island)

anno$chr <- gsub("chr", "", anno$chr) # remove chr from chromosome names
anno$chr <- gsub("X", "23", anno$chr) # replace chr X and Y with 23 and 24
anno$chr <- gsub("Y", "24", anno$chr)
anno$chr <- as.numeric(anno$chr) # set chr as numeric 

anno$UCSC_RefGene_Name	 <- as.factor(anno$UCSC_RefGene_Name) # set genes and strands as factors 
anno$strand <- as.factor(anno$strand)
anno[which(anno$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA #set missing gene names to NA 

anno$CpGs <- anno$Name
anno <- anno %>% dplyr::select(-Name)

cpg_anno <- left_join(highpip, anno, by = "CpGs")
fwrite(cpg_anno, file = "/6month/6month_anno_highpipcpgs.csv")