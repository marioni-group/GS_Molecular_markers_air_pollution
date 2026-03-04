#MAJA_methylation_dataprep
#phenotypes = 8 pollutants + smoking
#covariates = age, sex, Batch, logbmi, logunits, simdrank

#load libraries
library(tidyverse)
library(data.table)
library(future.apply)
library(arrow)
library(limma)
library(QCEWAS)
library(missMethods)
library("minfi")

library("optparse")
library("stringr")
library("imputeTS")
library("data.table")
library("limma")
library("foreach")
library("doParallel")

#Load covariates file for methylation data to create correct ID subset & order
cov <- fread("cov18512_nosmoksc_impt.csv")

#Load QC probes file
probes <- read.table("cpgs_tokeep.txt", header=F)

#Prepare methylation Data
# Set up input & output locations
datadir <- ""
methdir <- ""
localdir <- ""

#___________________________________________________________________________________________________________
#1. Load Data by Chromosome
#2. Subset to IDs with pollution data & covariates - 18512
#4. Subset probes to those passing QC
#5. Order match
#6. Impute if any NA values
#7. Regress for age, sex, batch, bmi, units and simd rank
#8. Write out CpGs
#9. Take residuals and scale across each CpG 
#   Check order
#10. Save out & exit 
#11. Fuse methylation files
#12. Make sure order is saved
#13. ?Save out full object

# Set up paralellizing
cores <- detectCores()  
cl <- makeCluster(3, outfile=paste0(datadir, "parallel_test.txt")) 
registerDoParallel(cl)

# Iterate per chromosome
foreach(i=1:22, .packages = "data.table") %dopar% { 

print(paste0("Working on chromosome ",i))

# Import 
  meth <- readRDS(paste0(localdir, "GS20k_chr", i, "_mvals.rds")) #Each chromosome saved out separately

  # Subset 
  meth <- meth[,which(colnames(meth) %in% cov$Sample_Sentrix_ID)] # Subset to those with complete relevant data 
  meth <- meth[which(rownames(meth) %in% probes$V1),] # Subset to probes passing QC     

  # Match order of IDs in phenotype and methylation file 
  meth <- meth[,match(cov$Sample_Sentrix_ID, colnames(meth))]


# Mean impute - cannot have NAs in final file 
  meth <- apply(meth, 1, meanimpute)
  meth <- t(meth) 


   # Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~age + as.factor(sex) +  as.factor(Batch) + logbmi + logunits + simd2009v2rank, data=cov)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth) 
  meth <- meth[!is.infinite(rowSums(meth)),]
  rm(fit.resid)
  gc()

 # Write out CpGs 
  cpgs <- as.data.frame(row.names(meth))
  SSID <- as.data.frame(colnames(meth))
  names(cpgs)[1] <- "CpG"
  names(SSID)[1] <- "SSID"
  fwrite(cpgs, paste0(datadir, methdir, "/GS_chr", i, "_cpgs.txt"),row.names=F)

 # Scale each CpG (column on transposed df) & transpose back    
  meth <- scale(t(meth)) 
  meth <- t(meth)
  gc()

  # Save out residualised file scaled
  meth <- data.table(meth, keep.rownames = TRUE)
  fwrite(meth, paste0(datadir, methdir, "/GS_chr", i, "_resid_mvals_scaled.txt"),row.names=F)  

 # Remove methylation object and clean up environment 
  rm(meth)
  gc()
}   


# End parallel
stopCluster(cl)

#Export files to zarr format and extract CpGs and SSIDs for matching/checking

# Generate file_list
file_list <- list.files(path = paste0(datadir, methdir), pattern = "_resid_mvals_scaled.txt", full.names = TRUE) #get files

#Load additional required libraries
library(arrow) #enables data manipulation without requiring full RAM
library(reticulate) #enables use of python
reticulate::py_install("zarr")
zarr <- import("zarr")
reticulate::py_install("pyarrow")
py_module_available("zarr")    # Should return TRUE
py_module_available("pyarrow") # Should return TRUE
import_builtins()

#Create function to process files in chunks and convert to zarr
process_files_to_zarr <- function(file_list, zarr_file) {
  combined_colnames <- NULL #to store column names
  row_mappings <- list() #to store row mappings (how each file relates to the rows)
  all_rownames <- list() #to store row names (CpGs)

  # Read the first file to get column information
  dt <- fread(file_list[1], data.table = TRUE)
  combined_colnames <- colnames(dt) #column names are SSIDs - should be the same for all files
  
  # Assuming the first column contains row names
  initial_rows <- 1e4
  zarr_store <- zarr$open_array(store = zarr$DirectoryStore(zarr_file),
                                mode = "w",
                                shape = tuple(as.integer(initial_rows), as.integer(ncol(dt) - 1)),
                                chunks = tuple(as.integer(1e3), as.integer(ncol(dt) - 1)),
                                dtype = "float32")

  row_start <- 0

  for (file in file_list) {
    dt <- fread(file, data.table = TRUE)
    
    if (!identical(combined_colnames, colnames(dt))) {
      stop("Column names do not match across files.")
    }
    
    # Extract and store row names
    rownames <- dt[[1]]
    all_rownames[[file]] <- rownames
    dt[[1]] <- NULL  # Remove the row names column for storage
    
    num_rows <- nrow(dt)
    row_end <- row_start + num_rows
    
    if (row_end > py_to_r(zarr_store$shape[1])) {
      zarr_store$resize(tuple(as.integer(row_end), as.integer(ncol(dt))))
    }
    
    # Use slices for proper indexing
    zarr_store[py_slice(as.integer(row_start), as.integer(row_end)),] <- as.matrix(dt)
    row_mappings[[file]] <- list(start = row_start + 1, end = row_end)
    row_start <- row_end
  }
  
  return(list(
    rownames_mapping = row_mappings,
    colnames = combined_colnames[-1],  # Exclude the row names column
    all_rownames = all_rownames  # Keep row names separately
  ))
}

file_list <- list.files(path = paste0(datadir, methdir), pattern = "_resid_mvals_scaled.txt", full.names = TRUE)
zarr_file <- paste0(datadir, methdir, "full_residsc_meth18512.zarr")
py_slice <- import_builtins()$slice
result <- process_files_to_zarr(file_list, zarr_file)
print(result$rownames_mapping)

#query rownames & colnames to check output
head(result$all_rownames)
head(result$colnames)


#Compare SSID orders to covariates file
cov <- fread("cov18512_nosmoksc_impt.csv")
identical(cov$Sample_Sentrix_ID, result$colnames) #TRUE
ssid_zarr_order <- as.data.frame(result$colnames)

#compare cpg output to combined cpg input

#combined cpg input
cpg_files <- list.files(path = paste0(datadir, methdir), pattern = "_cpgs.txt", full.names = TRUE)
cpg_order_check <- rbindlist(lapply(cpg_files, fread))


#cpg output
library(purrr)
class(result$all_rownames) 
cpgs_zarr <- as.data.frame(unlist(result$all_rownames, use.names = FALSE))
length(cpgs_zarr)  
identical(cpg_order_check$CpG, cpgs_zarr) 

#re-order covariates file if necessary -> not required

#save out orders
fwrite(ssid_zarr_order, file = "SSID_order_mvalinput.csv")
fwrite(cpgs_zarr, file = "cpg_order_mvalinput.csv")

#check shape of zarr file: want the shape to be rows = SSIDs, columns = CpGs
# Open the Zarr dataset in read-only mode

# import a Python module
zarr <- import("zarr")
z <- zarr$open("full_residsc_meth18512.zarr", mode="r")
shape <- z$shape
print(shape)


#Need to transpose the zarr file for use in MAJA - use python
library(reticulate)
reticulate::py_run_string("
import zarr
import numpy as np
zarr_data = zarr.open('full_residsc_meth18512.zarr', mode='r')

# Check the shape of the dataset
print('Shape of the Zarr dataset:', zarr_data.shape)

# Transpose and check shape
transposed_data = zarr_data[:].T  # This loads and transposes in memory
print('Transposed shape:', transposed_data.shape)
np.save('full_residsc_meth18512.npy', transposed_data)
z = zarr.array(transposed_data, store = 'full_residsc_meth18512_T.zarr', order = 'C')


# Open the Zarr dataset in read-only mode
zarr_data = zarr.open('full_residsc_meth18512_T.zarr', mode='r')

# Check the shape of the dataset
print('Shape of the Zarr dataset:', zarr_data.shape)
")

z <- zarr$open("full_residsc_meth18512_T.zarr", mode="r")
shape <- z$shape
print(shape) #check transposition has worked





check_cpgs <- fread("GS_chr1_cpgs.txt")
#Cpgs in column titled CpG
check_ids1 <- fread("GS_chr1_resid_mvals_scaled.txt", nrows = 3)
check_ids2 <- fread("GS_chr2_resid_mvals_scaled.txt", nrows = 3)
#CpGs are stored in first column of data

identical(colnames(check_ids1), colnames(check_ids2)) 
