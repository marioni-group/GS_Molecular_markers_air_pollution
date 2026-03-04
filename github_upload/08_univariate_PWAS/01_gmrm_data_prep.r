#Set up GMRMomi for all pollutants and smoking

#1. Create .bin file from pre-prepped protein data
#2. Convert phenotype to appropriate file type & save as separate files
#3. Create group index file
#4. Create .dim file (N = number individuals, M = number of markers)

#1. #GMRMomi with full pollutants: protein data prep
#________________________________________________________________________________________________________________________________

#Using protein data prepped for MAJA:
data <- fread("proteins_133prepped.csv")

#extract Proteins
proteins <- colnames(data)
proteins <- proteins[-1]
proteins <- as.data.frame(proteins)
fwrite(proteins, file = "/protein_order.csv")

#extract SSIDs
ids <- data$id
ids <- as.data.frame(ids)
fwrite(ids, file ="id_order.csv")

#remove from file and save out as bin file
#create .bin file

data <- data[,!"id"]
data[1:5, 1:5]
data <- unlist(data, use.names = FALSE)
data <- as.vector(data) 
head(data) #data should start with the first column as a vector

new = file("full_residsc_prot15314.bin", "wb")
writeBin(data, new)
close(new)

#check bin file with below
help(readBin)
check = file("full_residsc_prot15314.bin", "rb")
readBin(check, numeric(), n = 6)
close(check)



#2. Separate phenotype files
#________________________________________________________________________________________________________________________________
#Load pollution data - start with originals: all pollutants scaled, smoking, scaled

poll <- fread("/poll15314_sc.csv")

#check order of ids
identical(poll$id, ids$id) #TRUE

#add in FID/IID columns
##set up in required format
poll$FID <- poll$id
poll$IID <- poll$id
poll <- poll %>% select(-id)

dim(poll)
poll <- poll %>% dplyr::select(c(FID, IID, everything()))

location <- "/input_1yr/"
for(i in 3:11) { 
name <- as.character(names(poll)[i])
name_data <- poll %>%
  select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,".phen"), sep=" ", col.names = FALSE)
}

#check
check <- fread("input_1yr/PM25.phen")

#_______________________________________________________________________________
#Do the same for 6 month and 7yr pollution data
#6month
poll6mnth <- fread("/6month_avg_exp_foranalysis_ntr.csv")

#filter to ids for protein analysis
poll6mnth <- poll6mnth %>% filter(ID %in% poll$id)

#match order
poll6mnth <- poll6mnth[match(poll$id, poll6mnth$ID,)]
head(poll6mnth)
poll6mnth$ID <- as.character(poll6mnth$ID)
poll$id <- as.character(poll$id)
identical(poll6mnth$ID, poll$id) #TRUE

#scale
poll6mnth <- poll6mnth %>% mutate(across(where(is.numeric), scale))

#add in FID/IID columns
##set up in required format
poll6mnth$FID <- poll6mnth$ID
poll6mnth$IID <- poll6mnth$ID
poll6mnth <- poll6mnth %>% dplyr::select(-ID)

dim(poll6mnth)
poll6mnth <- poll6mnth %>% dplyr::select(c(FID, IID, everything()))

location <- ""
for(i in 3:10) { 
name <- as.character(names(poll6mnth)[i])
name_data <- poll6mnth %>%
  dplyr::select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,".phen"), sep=" ", col.names = FALSE)
}


#7 year data
poll7yr <- fread("/sevenyr_foranalysis_notsc.csv")
poll7yr <- poll7yr %>%
  dplyr::select(-c(raster_entry, standard_dev, year_range)) %>%
  pivot_wider(
    names_from = Pollutant,
    values_from = avg_7yrs
  )


#filter to ids for protein analysis
poll7yr <- poll7yr %>% filter(ID %in% poll$id)

#match order
poll7yr <- poll7yr[match(poll$id, poll7yr$ID),]
head(poll7yr)
poll7yr$ID <- as.character(poll7yr$ID)
poll$id <- as.character(poll$id)
identical(poll7yr$ID, poll$id) #TRUE

#scale
poll7yr <- poll7yr %>% mutate(across(where(is.numeric), scale))

#add in FID/IID columns
##set up in required format
poll7yr$FID <- poll7yr$ID
poll7yr$IID <- poll7yr$ID
poll7yr <- poll7yr %>% dplyr::select(-ID)

dim(poll7yr)
poll7yr <- poll7yr %>% dplyr::select(c(FID, IID, everything()))

#save out
location <- "input_7yr/"
for(i in 3:10) { 
name <- as.character(names(poll7yr)[i])
name_data <- poll7yr %>%
  dplyr::select(all_of(c("FID", "IID", name)))
fwrite(name_data, file = paste0(location, name,".phen"), sep=" ", col.names = FALSE)
}


#3. gri file
#________________________________________________________________________________________________________________________________
###Preparation of groups file
# 2 columns, first are probes, second is group assignation (cpgs vs. SNPs vs....)
# column length should = length of no. probes

example_group <- fread("example.gri")

#take protein labels 
head(proteins)
# Cpgs <- cpgs[-1]
# Cpgs <- as.data.frame(cpgs)

proteins$V2 <- c(1:133)
proteins$V3 <- c(0)
proteins <- proteins %>% select(-proteins)
fwrite(proteins, "full_gri.gri", sep = " ", col.names = FALSE)

#4. Dim file
#_______________________________________________________________________________________________________________________________
###Dim file
#x and y dimensions of original methylation matrix
#save with no column names, space separated

#example
example_dim <- fread("example.dim")
head(example_dim)


test_dim <- example_dim %>% mutate(V1 = 15314, V2 = 133)
fwrite(test_dim, "full_dim.dim", sep = " ", col.names = FALSE)


#5. Group mixture file
#________________________________________________________________________________________________________________________________
###Preparation of group mixture file

#Look at example
example_grm <- fread("grm.grm")
head(example_grm)


#Save out
fwrite(example_grm, file = "input_1yr/grm.grm", sep=" ", col.names = FALSE)
