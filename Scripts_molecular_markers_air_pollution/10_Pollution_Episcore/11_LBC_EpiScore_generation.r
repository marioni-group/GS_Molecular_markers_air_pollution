#Calculate Episcore in STRADL
library(tidyverse)
library(data.table)
library(utils)
library(here)
library(missMethods)
library(lumi)
#Load episcore weights
weights <- fread("cpg_weights.csv") 



#Load methylation data 
dat <- readRDS("LBC_betas_3489_bloodonly.rds")
wv1 <- fread("wv1_target36.csv")

#filter to wave1 data
dat[1:5, 1:5]
dim(dat) 

#Subset to target file
dat <- dat[,colnames(dat) %in% wv1$Basename]
dim(dat)

#convert to M values
dat2 <- beta2m(dat) #with CpGs as rows
dat2 <- t(dat2) #now CpGs are columns
dat <- dat2


#Transpose, scale, transpose back
dat <- scale(dat)

#Impute missing values
 sum(is.na(dat))


#Mean impute
dat <- as.data.frame(dat)
dat <- impute_mean(dat)
sum(is.na(dat))
dat <- t(dat)
# [1] 0



#create pollutant pollutant list
pollutants <- colnames(weights)[2:10]

#project episcores into wv1

#Calculate EpiScores_______________________________________________

#scale methylation data

#extract ssids
ssids <- colnames(dat)

# #Filter weights cpgs to those in data
weights <- weights[CpG %in% rownames(dat),]
dim(weights)

#Filter data to CpGs in weights
dat2 <- dat[rownames(dat) %in% weights$CpG,]
dim(dat2)


#match orders of cpg weights and methylation data
identical(weights$CpG, rownames(dat2)) 
all(weights$CpG %in% rownames(dat2)) 
all(rownames(dat2) %in% weights$CpG) 

weights <- weights[match(rownames(dat2), weights$CpG),]
identical(weights$CpG, rownames(dat2))
cpg_order <- rownames(dat2)


#Summarise CpGs in overlap between weights and 450k data

epi_lbccpg <- weights %>%
  pivot_longer(!CpG,
  names_to = "Trait",
  values_to = "Betas") %>%
  filter(!Betas == 0) %>%
  group_by(Trait) %>%
  summarise(
    weighted_cpgs = n(),
   positive_weights = sum(Betas > 0),
    negative_weights = sum(Betas < 0)) %>%
  ungroup()
fwrite(epi_lbccpg, file = "lbc_episcore_retained_cpgssummary.csv")



#Function to calculate episcore (multiplies CpGs by CpG weights)
calculate_episcore <- function (x,y) { 
  sum(x*y)
}

#Loop through all pollutants
#Methylation data should have rows as cpgs, SSIDs as columns

#loop requires data frame as input
dat2 <- as.data.frame(dat2)

episcore_results <- list()  #initialise empty list


for (pollutant in pollutants) {
 y <- as.data.frame(weights[[pollutant]]) #extracts CpG weights for protein 
 print(pollutant)
episcore_results[[pollutant]] <- purrr::map(dat2, ~(calculate_episcore(.x,y)),
                     .id = "SSID") #
}

epi_df <- episcore_results %>% lmap(bind_rows, .id = "pollutant") %>% bind_rows()  
head(epi_df)
fwrite(epi_df, file = "lbc_wv1_episcores.csv")