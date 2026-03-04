#Compare EpiScore and pollution exposure

library(tidyverse)
library(data.table)
library(performance)

wv2 <- fread("stradl_wv2_episcores_unrelated.csv")
wv3 <- fread("stradl_wv3_episcores_unrelated.csv") 
pollution <- fread("stradl_365_av_sc.csv")

wv2 <- dcast(melt(wv2, id.vars = "pollutant"), variable ~ pollutant)
wv3 <- dcast(melt(wv3, id.vars = "pollutant"), variable ~ pollutant)
episcore <- rbind(wv2, wv3)
head(episcore)
episcore <- setDT(episcore)
episcore$variable <- as.character(episcore$variable)

#Match IDs and orders

#load id match file
wv2_ids <- readRDS("w2-stradl-pData.rds")
wv2_ids <- wv2_ids %>% select(Sample_Name, Sample_Sentrix_ID)
wv2_ids$Sample_Name <- gsub("S", "", wv2_ids$Sample_Name)

wv3_ids <- fread("pdata-with-ST_age.csv")
wv3_ids <- wv3_ids %>% select(Sample_Name, Sample_Sentrix_ID)
ids <- rbind(wv2_ids, wv3_ids)
ids <- ids %>% rename(ID = Sample_Name)
pollution$ID <- as.character(pollution$ID)
pollution2 <- left_join(pollution, ids, by = "ID")

pollution2 <- pollution2 %>% 
  drop_na() 

#
pollution3 <- pollution2 %>%
  filter(Sample_Sentrix_ID %in% episcore$variable) 


#scale episcore
str(episcore)
episcore <- episcore %>% filter(variable %in% pollution3$Sample_Sentrix_ID)
episcore <- episcore %>% mutate(across(where(is.numeric), scale))

#match orders of pollution & episcore
episcore <- episcore[match(pollution3$Sample_Sentrix_ID, episcore$variable),]
identical(episcore$variable, pollution3$Sample_Sentrix_ID)

#Prep for plotting/regression
episcore <- episcore %>% mutate(type = "episcore") %>% rename(Sample_Sentrix_ID = variable) %>% select(-logsmok)
episcore <- episcore %>% pivot_longer(.,
    !c(Sample_Sentrix_ID, type),
    names_to = "pollutant",
    values_to = "value")

measured <- pollution3 %>% mutate(type = "modelled") %>% select(-ID)
measured <- measured %>% pivot_longer(.,
     !c(Sample_Sentrix_ID, type),
     names_to = "pollutant",
     values_to = "value")

#combine
plot_df <- rbind(episcore, measured)
plot_df1 <- plot_df %>% pivot_wider(.,
                                    names_from = "type",
                                   values_from = "value")
head(plot_df1)                                   
sum(is.na(plot_df1))

#Function calculating r
cor_predictr <-   function(x){
  return(
          broom::tidy(cor.test(x$episcore, x$modelled, method = "pearson"))) 
}

#on 693 individuals
correlation_results <- plot_df1 %>%
  group_by(pollutant) %>%
  do(cor_predictr(.))
correlation_results <- correlation_results %>% arrange(desc(estimate))


fwrite(correlation_results, file = "cor_stradlunrelated_episcore.csv")


