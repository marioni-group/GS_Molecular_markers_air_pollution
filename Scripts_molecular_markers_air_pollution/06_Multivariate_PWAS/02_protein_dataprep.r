#2. Regress covaraites out of protein data
library(tidyverse)
library(data.table)
library(here)


#Load data
proteins <- readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")
cov <- fread("cov15314_nosmoksc_impt.csv")
pollution <- fread("annual_avg_exp_foranalysis_scaled.csv")

#check ids and match orders
dim(cov) 

#filter
proteins2 <- proteins %>% 
  dplyr::filter(id %in% cov$id)

#match
proteins2 <- proteins2[match(cov$id, proteins2$id),]
identical(proteins2$id, cov$id) #TRUE

#Remove protein groups
head(proteins2)

proteins3 <- proteins2[, !grepl("\\.", colnames(proteins2))]



#Take residuals of protein data for covariates
#make protein list
protein_list <- colnames(proteins3)
protein_list <-protein_list[-1]

#take residuals
residualised = proteins3
for (protein in protein_list) {
  print(protein)
  residualised[[protein]] <- paste(resid(lm(residualised[[protein]] ~ cov$age + as.factor(cov$sex) + cov$logbmi + cov$logunits + cov$simd2009v2rank, 
   data = residualised, 
   na.action = na.exclude)))
}

#residualised data frame contains columns as type character -> need to be converted
residualised[1:5, 1:5]
residualised_2 <- residualised
residualised_2[2:134] <- lapply(residualised[2:134], as.numeric)

residualised_scaled_2 = residualised_2 
for (protein in protein_list) {
    print(protein)
    residualised_scaled_2[[protein]] <- paste(scale(residualised_2[[protein]]))
}



residualised_scaled_2[1:5, 1:5]
fwrite(residualised_scaled_2, file = "proteins_133prepped.csv")


#check order with pollution and covariates input
cov <- fread("cov15314_sc_impt.csv")
identical(cov$id, residualised_scaled_2$id) #TRUE


#convert to zarr 
#Load additional required libraries
library(arrow) #enables data manipulation without requiring full RAM
library(reticulate) #enables use of python
reticulate::py_install("zarr")
zarr <- import("zarr")
# numcodecs <- import("numcodecs")
reticulate::py_install("pyarrow")
py_module_available("zarr")    # Should return TRUE
py_module_available("pyarrow") # Should return TRUE
import_builtins()
pyreadr <- import("pyreadr")


data <- residualised_scaled_2 %>% dplyr::select(-id)


data[1:5, 1:5]
str(data)
data <- purrr::map_df(data, as.numeric) 

str(data)
matrix_data <- as.matrix(data)
chunk_shape <- c(10L,10L)
store <- zarr$array(matrix_data, 
                    store = 'protein_covcorr_0401.zarr', 
                    chunks = tuple(chunk_shape),
                    config = list(order = "C"),
                    )
zarr_data = zarr$open('protein_covcorr_0401.zarr', mode='r')
zarr_data$shape



