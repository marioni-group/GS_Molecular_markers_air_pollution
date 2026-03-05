#Gene ontology/enrichment of significant CpGs from MAJA pollution ewas
library(minfi)
library(data.table)
library(tidyverse)
library(missMethyl)



maja_results <- fread("majaewas_highpip_anno_0812_newthresh.csv")
maja_results <- maja_results %>% filter(pip > 0.95 & association2 == "yes")
included_probes <- fread("cpgs_tokeep.txt", header = FALSE) 
epic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

cpgs <- maja_results$CpG
bg <- included_probes$V1

enrich <- gometh(cpgs,
      bg,
      collection = c("GO", "KEGG"),
      array.type = "EPIC",
      plot.bias = TRUE,
      prior.prob = TRUE,
      anno = epic,
      equiv.cpg = TRUE,
      fract.counts = TRUE,
      genomic.features = "ALL",
      sig.genes = TRUE
)

head(enrich)
enrich %>% dplyr::filter(FDR < 0.05)
limma::topGO(enrich)
fwrite(enrich, file = "enrichment_gometh_0912_newthres.csv")

