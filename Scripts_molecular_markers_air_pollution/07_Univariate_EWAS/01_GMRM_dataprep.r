#Set up GMRMomi for all pollutants and smoking

#1. Convert prepped methylation data to .bin
#2. Convert phenotype to appropriate file type & save as separate files
#3. Create group index file
#4. Create .dim file (N = number individuals, M = number of markers)

#1. #GMRMomi with full pollutants: methylation prep
#________________________________________________________________________________________________________________________________

#From original files
datadir <- ""
methdir <- ""

# Generate file_list
file_list <- list.files(path = paste0(datadir, methdir), pattern = "_resid_mvals_scaled.txt", full.names = TRUE) #get files
file_list <- file_list[order(file_list)]
data <- rbindlist(lapply(file_list,fread)) #fuse methylation data
gc()

# Export fused methylation file - scaled 
data <- data.table(data, keep.rownames = TRUE)
data[1:5, 1:5]
#tranpose so IDs = rows, cpgs = columns
data <- data.table::transpose(data, keep.names = "rn", make.names = "rn") 
gc()

#extract CpGs
cpgs <- colnames(data)
cpgs <- cpgs[-1]
cpgs <- as.data.frame(cpgs)
fwrite(cpgs, file = "cpg_order_methinput.csv")

#extract SSIDs
ssids <- data$rn
ssids <- as.data.frame(ssids)
fwrite(ssids, file ="ssid_order_methinput.csv")

#remove from file and save out as bin file
#create .bin file
# data <- data %>% select(-rn) #meth should be in the format of IDs(rows) x CpGs(columns)
data <- data[,!"rn"]
data[1:5, 1:5]
data <- unlist(data, use.names = FALSE)
data <- as.vector(data) 
head(data) #data should start with the first column as a vector

new = file("full_residsc_meth18512_T2.bin", "wb")
writeBin(data, new)
close(new)

#check bin file with below
help(readBin)
check = file("full_residsc_meth18512_T2.bin", "rb")
readBin(check, numeric(), n = 6)
close(check)



#2. Separate phenotype files
#________________________________________________________________________________________________________________________________
#Load pollution data - start with originals: all pollutants scaled, smoking, scaled

poll <- fread("full_poll_ids_newphenorder_2905.csv")

#add in FID/IID columns
test <- fread("input/NO.phen")

##set up in required format
poll$FID <- poll$id
poll$IID <- poll$id
poll <- poll %>% select(-id)

dim(poll)
poll <- poll %>% dplyr::select(c(FID, IID, everything()))
poll <- poll %>% select(-Sample_Sentrix_ID)
location <- ""
for(i in 3:11) { 
name <- as.character(names(poll)[i])
name_data <- poll %>%
  select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,".phen"), sep=" ", col.names = FALSE)
}

#3. gri file
#________________________________________________________________________________________________________________________________
###Preparation of groups file
# 2 columns, first are probes, second is group assignation (cpgs vs. SNPs vs....)
# column length should = length of no. probes

example_group <- fread("example.gri")
cpgs <- fread("cpg_order_methinput.csv")
#take CpG labels from prepped methylation data
head(cpgs)
# Cpgs <- cpgs[-1]
Cpgs <- as.data.frame(cpgs)
Cpgs$V2 <- c(0)


fwrite(Cpgs, "full_gri.gri", sep = " ", col.names = FALSE)
fwrite(Cpgs, "full_gri.gri", sep=" ", col.names = FALSE)
#4. Dim file
#_______________________________________________________________________________________________________________________________
###Dim file
#x and y dimensions of original methylation matrix
#save with no column names, space separated

#example
example_dim <- fread("example.dim")
head(example_dim)

#check numbers!!
test_dim <- example_dim %>% mutate(V1 = 18512, V2 = 752722)
fwrite(test_dim, "full_dim.dim", sep = " ", col.names = FALSE)


#5. Group mixture file
#________________________________________________________________________________________________________________________________
###Preparation of group mixture file

#Look at example
example_grm <- fread("grm.grm")
head(example_grm)


#Save out
fwrite(example_grm, file = "grm.grm", sep="\t", col.names = FALSE)

#_______________________________________________________________________________________________________
#Also set up for 6month and 7yr data

#GMRM set up for 7-year & 6 month data

library(tidyverse)
library(data.table)

#6 month
poll <- fread("6month_avg_exp_foranalysis.csv")

#filter to N for ewas
comp <- fread("full_poll_ids_newphenorder_2905.csv")

poll <- poll %>% filter(ID %in% comp$id)
identical(poll$ID, comp$id)

poll <- poll[match(comp$id, poll$ID),]
identical(poll$ID, comp$id) #TRUE

##set up in required format
poll$FID <- poll$ID
poll$IID <- poll$ID
poll <- poll %>% select(-ID)

dim(poll)
poll <- poll %>% dplyr::select(c(FID, IID, everything()))
poll$FID <- as.character(poll$FID) 
poll$IID <- as.character(poll$IID)
poll <- poll %>% mutate(across(where(is.numeric), scale))


location <- ""
for(i in 3:10) { 
name <- as.character(names(poll)[i])
name_data <- poll %>%
  select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,"6month.phen"), sep=" ", col.names = FALSE)
}

#7 years
poll <- fread("sevenyr_foranalysis.csv")
poll <- poll %>%
  select(-c(raster_entry, standard_dev, year_range)) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = avg_7yrs
  )

poll <- poll %>% filter(ID %in% comp$id)
identical(poll$ID, comp$id)

poll <- poll[match(comp$id, poll$ID),]
identical(poll$ID, comp$id) #TRUE

##set up in required format
poll$FID <- poll$ID
poll$IID <- poll$ID
poll <- poll %>% select(-ID)

dim(poll)
poll <- poll %>% dplyr::select(c(FID, IID, everything()))
poll$FID <- as.character(poll$FID) 
poll$IID <- as.character(poll$IID)
poll <- poll %>% mutate(across(where(is.numeric), scale))

location <- ""
for(i in 3:10) { 
name <- as.character(names(poll)[i])
name_data <- poll %>%
  select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,"7yr.phen"), sep=" ", col.names = FALSE)
}
