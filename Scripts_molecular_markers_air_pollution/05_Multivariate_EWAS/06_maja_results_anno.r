#Annotate Maja pollution ewas results

library(tidyverse)
library(data.table)
library(qqman)
library(UpSetR)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


results <- fread("maja_fullresults_1307_newthres.csv")
results_hp <- results %>% dplyr::filter(pip > 0.95)

#Summarise non-zero betas
nonzerocpgs <- results %>%
  group_by(Trait) %>%
  summarise(Count = sum(Betas != 0))
nonzero <- results %>%
  filter(Betas !=0 & pip > 0)  

epic_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- data.frame(epic_anno)
setDT(anno, keep.rownames = "CpG")
colnames(anno)
anno2 <- anno %>% dplyr::select(CpG, chr, pos, strand,Relation_to_Island,UCSC_RefGene_Name,Methyl450_Loci)

#For easy inspection
results_hpanno <- left_join(results_hp, anno2, by = "CpG")

fwrite(results_hpanno, file = "/maja_highpip_anno_0725")


results_hpanno2 <- results_hpanno %>%
 group_by(CpG) %>%
 dplyr::filter(any(association == "yes")) %>%
 ungroup() #all remain
fwrite(results_hpanno2, file = "/majaewas_highpip_anno_0812_newthresh.csv")

#for text
hp <- results_hpanno %>% filter(association2 == "yes")
dim(hp)
head(hp)
unique(hp$UCSC_RefGene_Name)
unique(hp$CpG)

hp %>% filter(Trait == "logsmok") %>% select(CpG) %>% unique()


#For plotting
results_anno <- left_join(results, anno2, by = "CpG")
length(unique(results_anno$CpG))

results_hpanno_upset <- results_anno %>%
  dplyr::filter(pip > 0.95 & association == "yes")
unique(results_hpanno_upset$UCSC_RefGene_Name)




        