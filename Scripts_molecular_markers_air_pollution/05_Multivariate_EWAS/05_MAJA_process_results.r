#Maja protein-pollution process results

library(tidyverse)
library(data.table)
library(utils)
library(here)

betas<- fread(cmd = "unzip -cqp /mean_beta.csv.zip")
betas <- betas[-1,]
prob <- fread("mean_prob.txt")
prob %>% dplyr::filter(V1 > 0.95) %>% dim() 

#add in proteins
cpgs <- fread("cpg_order_mvalinput.csv")
colnames(cpgs) <- "CpG"
probs <- cbind(cpgs, prob)
probs <- probs %>% dplyr::rename(pip = V1)
highpip <- probs %>% dplyr::filter(pip > 0.95) 

#add proteins to betas
betas <- cbind(cpgs, betas)


#variance
variance <- fread(cmd = "unzip -cqp /var_beta.csv.zip")
head(variance)
variance <- variance[-1,]
variance <- cbind(cpgs, variance)

#Create results for all input/output
results_full <- left_join(betas, probs, by = c("CpG"))
results_full <- results_full %>%
  dplyr::select(CpG, pip, everything())

#Extract pollutant names
poll <- fread("full_poll_ids_newphenorder_2905.csv")
poll <- poll %>% dplyr::select(-c(id, Sample_Sentrix_ID))
trait_list <- colnames(poll)[1:9]

#Rename results output
colnames(results_full)[3:11] <- trait_list 

results_full2 <- results_full %>%
  pivot_longer(!c(CpG, pip),
               names_to = "Trait",
               values_to = "Betas")

#Add pollutants to variance output
colnames(variance)[2:10] <- trait_list
variance <- variance %>%
  pivot_longer(
    !CpG,
    names_to = "Trait",
    values_to = "Variance"
  )

#Combine variance with other results
results_full2 <- left_join(results_full2, variance, by = c("CpG", "Trait"))

#Calculate standard deviations
results_full2 <- results_full2 %>%
    mutate(stdev = sqrt(Variance))

#Identify associations based on 1SD and 2SD credible intervals
results_full2 <- results_full2 %>%
    mutate(lower = Betas - stdev,
        upper = Betas + stdev,
        association = ifelse(pip > 0.95 & lower < 0 & upper < 0, "yes",
                        ifelse(pip > 0.95 & lower > 0 & upper > 0, "yes", "no")),
        associated_direction = ifelse(Betas < 0 & association == "yes", "Negative",
                                        ifelse(Betas > 0 & association == "yes", "Positive", NA)))


#ADD in 2SD cut off for results
results_full2 <- fread("maja_fullresults_1307.csv")
results_full3 <- results_full2 %>%
  mutate(lower2 = Betas - 1.96*stdev,
          upper2 = Betas + 1.96*stdev,
                  association2 = ifelse(pip > 0.95 & lower2 < 0 & upper2 < 0, "yes",
                        ifelse(pip > 0.95 & lower2 > 0 & upper2 > 0, "yes", "no")))

#Explore associations
table(results_full3$association)
table(results_full3$association2)
results_full3 %>% filter(association2 == "yes") 
results_full3 %>% filter(association2 == "yes" & Trait == "logsmok") 
results_full3 %>% filter(association2 == "yes" & !Trait == "logsmok") 
fwrite(results_full3, file = "maja_fullresults_1307_newthres.csv")
results_full3 <- fread("maja_fullresults_1307_newthres.csv")

#________________________________________________________________________________________
#Association counts for CI of 1SD ***check supplementary***
highpip_results <- results_full2 %>% 
  dplyr::filter(pip > 0.95) %>%
  dplyr::select(CpG, pip, Trait, Betas, association, associated_direction)

table(highpip_results$association, highpip_results$Trait)

fwrite(highpip_results, file = "maja_highpipresults_1307.csv")
highpip_results <- fread("maja_highpipresults_1307.csv")


highpip_plot <- highpip_results %>%
  dplyr::filter(association == "yes")

trait_counts <- highpip_plot %>%
  dplyr::count(Trait, sort = TRUE)

trait_counts %>% 
  ggplot() +
  geom_col(aes(x = Trait, y = n))
trait_counts$Trait <- forcats::fct_reorder(trait_counts$Trait, desc(trait_counts$n))  

p1 <- trait_counts %>%
  ggplot(aes(x = Trait, y = n)) +
  geom_col(fill = "#85A6F8", colour = "black") +
  xlab("Pollutant") +
  ylab("Number associated proteins") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5)),
        legend.position="none") +
  coord_cartesian(expand = FALSE, xlim = c(0, NA), ylim = c(0, 25))
ggsave(p1, file = "assoc_bypoll_0725.pdf")



summary_df <- highpip_results %>%
  filter(association == "yes") %>%
  group_by(CpG) %>%
  summarise(associated_traits = paste(Trait, collapse = ","),
            trait_count = n()) %>%
  arrange(desc(trait_count))

fwrite(summary_df, file = "traits_by_cpg_0925.csv")


CpG_counts$CpG <- forcats::fct_reorder(CpG_counts$CpG, desc(CpG_counts$n)) 

p2 <- CpG_counts %>%
  ggplot(aes(x = CpG, y = n)) +
  geom_col(fill = "#85A6F8", colour = "black") +
  xlab("CpG") +
  ylab("Number associated pollutants") +
  scale_y_continuous(breaks = c(0:15)) +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.x = element_blank(),
        legend.position="none") +
  coord_cartesian(expand = FALSE, xlim = c(0, NA), ylim = c(0, 10))
ggsave(p2, file = "assoc_byCpG_0725.pdf" )

summary_df <- results_full2 %>%
  filter(association == "yes") %>%
  group_by(CpG) %>%
  summarise(associated_traits = paste(Trait, collapse = ","),
            trait_count = n())

summary_df <- summary_df %>% arrange(desc(trait_count))
fwrite(summary_df, file = "maja_summary_0725.csv")


#________________________________________________________________________________________
#Association counts for CI of 2SD

highpip_results2 <- results_full3 %>% 
  dplyr::filter(pip > 0.95) %>%
  dplyr::select(CpG, pip, Trait, Betas, association2, associated_direction)

table(highpip_results2$association2, highpip_results2$Trait)

fwrite(highpip_results2, file = "maja_highpipresults_1307_newthresh.csv")

#Counts for text
df <- highpip_results2 %>% filter(association2 == "yes")
dim(df)
length(unique(df$CpG))

highpip_plot <- highpip_results2 %>%
  dplyr::filter(association2 == "yes")

trait_counts <- highpip_plot %>%
  dplyr::count(Trait, sort = TRUE)

trait_counts %>% 
  ggplot() +
  geom_col(aes(x = Trait, y = n))
trait_counts$Trait <- forcats::fct_reorder(trait_counts$Trait, desc(trait_counts$n))  

p1 <- trait_counts %>%
  ggplot(aes(x = Trait, y = n)) +
  geom_col(fill = "#85A6F8", colour = "black") +
  xlab("Pollutant") +
  ylab("Number associated proteins") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5)),
        legend.position="none") +
  coord_cartesian(expand = FALSE, xlim = c(0, NA), ylim = c(0, 25))
ggsave(p1, file = "assoc_bypoll_0812_newthresh.pdf")
  

summary_df <- highpip_results2 %>% #CpGs & associated traits
  filter(association2 == "yes") %>%
  group_by(CpG) %>%
  summarise(associated_traits = paste(Trait, collapse = ","),
            trait_count = n()) %>%
  arrange(desc(trait_count))

fwrite(summary_df, file = "traits_by_cpg_0812_newthresh.csv")

CpG_counts <- highpip_plot %>%
  dplyr::count(CpG, sort = TRUE)


CpG_counts$CpG <- forcats::fct_reorder(CpG_counts$CpG, desc(CpG_counts$n)) 

p2 <- CpG_counts %>%
  ggplot(aes(x = CpG, y = n)) +
  geom_col(fill = "#85A6F8", colour = "black") +
  xlab("CpG") +
  ylab("Number associated pollutants") +
  scale_y_continuous(breaks = c(0:15)) +
  theme_bw() +
  theme(axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.x = element_blank(),
        legend.position="none") +
  coord_cartesian(expand = FALSE, xlim = c(0, NA), ylim = c(0, 10))
ggsave(p2, file = "assoc_byCpG_0812_newthresh.pdf" )


#Variance explained___________________________________________________________________

df <- fread("/mean_V.txt")

dim(df)
df2 <- fread("var_V.txt")

names(df) <- trait_list
ids <- names(df)
var <- diag(as.matrix(df))
var_var <- diag(as.matrix(df2))


#make df with variance and variance of variance 
var_df <- data.frame(Measure = ids, Variance = var)
var_var_df <- data.frame(Measure = ids, Variance_of_variance = var_var)

#merge together 
all <- full_join(var_df, var_var_df, by = "Measure")

#multiply variance and variance of variance by 100
all$Variance <- all$Variance * 100
all$Variance_of_variance <- all$Variance_of_variance * 100

#create upper and lower variance limit 
all$lower <- all$Variance - all$Variance_of_variance
all$upper <- all$Variance + all$Variance_of_variance

#Save out
fwrite(all, file = "ewas_variances_0725.csv")


#plot
all$Measure <-forcats::fct_reorder(all$Measure, dplyr::desc(all$Variance))

library(paletteer)
clist <- paletteer_d("MetBrewer::Tam")
clist <- c( "#910000FF", clist)
#reorder 
clist <- c(logsmok = "#910000FF",
                                 NO = "#FFD353FF",
                                 NO2 = "#FFB242FF",
                                 NO3_C = "#EF8737FF",
                                 NO3_F = "#DE4F33FF",
                                 O3 = "#BB292CFF",
                                 PM10 = "#9F2D55FF",
                                 PM25 = "#62205FFF",
                                 SO2 = "#341648FF")

all %>%
  ggplot(aes(x = Measure, y = Variance, fill = Measure)) +
  scale_fill_manual(values = clist) +
  geom_col() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  theme_bw() +
  labs(x = "Pollutant")
ggsave("ewas_variances_plot.png")

#_______________________________________________________
#Covariances
#create covar df by extracting the upper triangle of dataset (lower triangle also works)
covar <- as.matrix(df)
rownames(covar) <- colnames(covar)
covariances <- covar[upper.tri(covar)]

#make row names ans column names ids 
row_names <- rownames(covar)
col_names <- colnames(covar)

#create labels for covariances 
labels <- outer(row_names, col_names, FUN = paste, sep = "-")[upper.tri(covar)]

#create dataframe with covariances and labels 
labeled_covariances <- data.frame(Covariances = covariances, Labels = labels)

#variance of covariance matrix 
covar_var <- as.matrix(df2)

#set row names and column names to be the same 
rownames(covar_var) <- trait_list
colnames(covar_var) <- trait_list

#extract variance of covariance (upper triangle or lower triangle )
covariances_var <- covar_var[upper.tri(covar_var)]

#set row names and column names 
row_names <- rownames(covar_var)
col_names <- colnames(covar_var)

#create labels for variance of covariance matrix 
labels <- outer(row_names, col_names, FUN = paste, sep = "-")[upper.tri(covar_var)]

#make dataframe with labels and variances of covaiances 
labeled_covariances_var <- data.frame(Covariances_var = covariances_var, Labels = labels)

#merge together 
all_covar <- full_join(labeled_covariances, labeled_covariances_var, by = "Labels")

#multiple covariance and variance of covariance by 100
all_covar$Covariances <- all_covar$Covariances * 100 
all_covar$Covariances_var <- all_covar$Covariances_var * 100

min(abs(all_covar$Covariances))

max(abs(all_covar$Covariances))


#save out
fwrite(all_covar, file = "maja_covariances_0725.csv")


#Plot
all$Measure <- fct_reorder(all$Measure, desc(all$Variance))

v1  <- all %>%
  ggplot(aes(x = Measure, y = Variance, colour = Measure)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  theme_bw() +
  labs(x = "Pollutant")
ggsave(v1, file = "variance_0725.pdf", width = 20, height = 10, units = "cm")  
 
all_covar$Labels <- fct_reorder(all_covar$Labels, desc(all_covar$Covariances)) 
v2 <-  all_covar %>%
  ggplot(aes(x = Labels, y = Covariances, colour = Labels)) + 
  geom_point(size = 1) +
  theme_bw() + 
  labs(x = "Pollutants") +
  theme(legend.position = "none",
         axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
ggsave(v2, file = "covariance_0725.pdf", width = 20, height = 10, units = "cm")

#___________________________________________________________________________________


