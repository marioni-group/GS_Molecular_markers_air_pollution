#Calculate Episcore in STRADL using all methylation data

#Load episcore weights
weights <- fread("cpg_weights.csv") 

#Using methylation data with everyone (including related)
w2 <- readRDS("w2-unrelated-STRADLonly-MVALS.rds") 
w3 <- readRDS("w3-unrelated-STRADLonly-MVALS.rds") 


#create pollutant pollutant list
pollutants <- colnames(weights)[2:10]

#project episcores into wv2 

#Calculate EpiScores_______________________________________________

#scale methylation data
#transpose so cpgs are columns
#extract ssids
ssids <- colnames(w2)

#Filter cpgs to those in weights
w2 <- w2[rownames(w2) %in% weights$CpG,]
dim(w2)

#Transpose, scale, transpose back
w2 <- t(w2)
w2 <- scale(w2)
w2 <- t(w2)


#match orders of cpg weights and methylation data
identical(weights$CpG, rownames(w2)) #FALSE
weights <- weights[match(rownames(w2), weights$CpG),]
identical(weights$CpG, rownames(w2)) #TRUE
cpg_order <- rownames(w2)


#Function to calculate episcore (multiplies CpGs by CpG weights)
calculate_episcore <- function (x,y) { 
  sum(x*y)
}

#Loop through all proteins
#Methylation data should have rows as cpgs, SSIDs as columns

#loop requires data frame as input
w2 <- as.data.frame(w2)

episcore_results <- list()  #initialise empty list


for (pollutant in pollutants) {
 y <- as.data.frame(weights[[pollutant]]) #extracts CpG weights for protein 
 print(pollutant)
episcore_results[[pollutant]] <- purrr::map(w2, ~(calculate_episcore(.x,y)),
                     .id = "SSID") #
}

epi_df <- episcore_results %>% lmap(bind_rows, .id = "pollutant") %>% bind_rows()  
head(epi_df)
fwrite(epi_df, file = "stradl_wv2_episcores_unrelated.csv")


#project episcores into wv3 

#Calculate EpiScores_______________________________________________

#scale methylation data
#transpose so cpgs are columns
#extract ssids
ssids <- colnames(w3)

#Filter cpgs to those in weights
w3 <- w3[rownames(w3) %in% weights$CpG,]
dim(w3)

#Transpose, scale, transpose back
w3 <- t(w3)
w3 <- scale(w3)
w3 <- t(w3)


#match orders of cpg weights and methylation data
identical(weights$CpG, rownames(w3))
weights <- weights[match(rownames(w3), weights$CpG),]
identical(weights$CpG, rownames(w3))
cpg_order <- rownames(w3)


#Function to calculate episcore (multiplies CpGs by CpG weights)
calculate_episcore <- function (x,y) { 
  sum(x*y)
}

#Loop through all proteins
#Methylation data should have rows as cpgs, SSIDs as columns

#loop requires data frame as input
w3 <- as.data.frame(w3)

episcore_results <- list()  #initialise empty list


for (pollutant in pollutants) {
 y <- as.data.frame(weights[[pollutant]]) #extracts CpG weights for protein 
 print(pollutant)
episcore_results[[pollutant]] <- purrr::map(w3, ~(calculate_episcore(.x,y)),
                     .id = "SSID") #
}

epi_df <- episcore_results %>% lmap(bind_rows, .id = "pollutant") %>% bind_rows()  
head(epi_df)
fwrite(epi_df, file = "stradl_wv3_episcores_unrelated.csv")
