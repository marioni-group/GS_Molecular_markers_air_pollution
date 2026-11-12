install.packages(c("mediation", "tidyverse"))
library(mediation)
library(tidyverse)
library(coxme)
library(data.table)

#Data set-up 
# subset to CpG/protein overlap as for protein EWAS
# - extract M values for CpG 
# - extract protein data 
#set-up covariates 
    # MAJA takes into account white cell proportions and genetic structure, so need to adjust models for these
    #Covariates in model: age + as.factor(sex) +  as.factor(Batch) + logbmi + logunits + simd2009v2rank


#Path: (CpG as Mediator): Pollutant ->CpG Site -> Protein.

# # 1. Mediator Model (Exposure -> Mediator)
# model_mediator <- lm(cpg_site ~ pollutant + age + sex + smoke_status + cell proportions + kinship, data = my_data)

# # 2. Outcome Model (Exposure + Mediator -> Outcome)
# model_outcome <- lm(protein ~ pollutant + cpg_site + age + sex + smoke_status + cell proportions + kinship, data = my_data)


#1. Prep data: 
    #CpG data: take residuals of the model: lme(CpG ~ age + sex + smoking, alcohol, bmi, cell proportions, simd + (1|id)) -> standardise
    #Protein data: take residuals of the model: lme(protein ~ age + sex + smoking, alcohol, bmi, cell proportions, simd + (1|id)) -> standardise
    #Covariates: smoking, alcohol, bmi, cell proportions, simd -> transform where necessary & standardise (cell proportions for CpG model only?)

#CpG = cg26561082, chromosome 10
#protein = P19827

#1__________________________________________________________________
library(lme4)
library(VIM)


localdir <- ""

#Load Data
cg_dat <- readRDS(paste0(localdir, "Chromosomes/GS20k_chr10_mvals.rds")) #Each chromosome saved out separately 
wcc <- readRDS(paste0(localdir, "GS20k_Targets_18869.rds"))
prot <- readRDS("GS_Proteins.rds")
keep <- fread("2024-10-31_meth_ids.csv")

#Residualise cpg site 
    #1. Select CpG data in question
    #2. Filter to IDs in protein data 
    #3. Create covariates df -> filter to BMI as in EWAS 
    #4. Filter CpG data by filtered DF
    #5. Residualise CpG data & standardise residuals
    #6. Save out result

#1. 
cgdat <- t(cg_dat)
cgdat <- data.frame(cgdat)
head(row.names(cgdat))
head(colnames(cgdat))
cgdat <- cgdat %>% dplyr::select("cg26561082") %>%
   mutate(Sample_Sentrix_ID = row.names(cgdat))

#2.
ids <- wcc %>% dplyr::select(Sample_Name, Sample_Sentrix_ID) %>%
        dplyr::rename(id = Sample_Name)
cgdat <- left_join(cgdat, ids, by = "Sample_Sentrix_ID")        
cgdat <- cgdat %>% filter(id %in% prot$id) #14671
cgdat <- cgdat %>% filter(id %in% keep$V1) #14665

#3.
cov <- fread("cov18512_noimpt.csv")

#add wcc
wcc2 <- wcc %>% dplyr::select(Sample_Name, Bcell, CD4T, CD8T, Eos, Mono, Neu, NK) %>%
       dplyr::rename(id = Sample_Name)
cov$id <- as.character(cov$id)
cov <- dplyr::left_join(cov, wcc2, by = "id")

cov2 <- cov %>%
dplyr::filter(bmi >=16 & bmi <=50 | is.na(bmi))

#Filter data to get the minimal overlapping sample
ids <- intersect(cgdat$id, prot$id) 
cov2 <- cov2 %>% filter(id %in% ids) 
cgdat <- cgdat %>% filter(id %in% cov2$id)
prot <- prot %>% filter(id %in% cov2$id)



#log transform
cov2 <- cov2 %>%
  mutate(logbmi = log(bmi),
         logunits = log(units + 1),
         logsmok = log(pack_years + 1))
sapply(cov2, function(y) sum(length(which(is.na(y)))))


#6. Impute missing data (KNN)
#Prep for scaling and knn imputation
 #Female is 1
cov2 <- cov2 %>%
  mutate(sex = ifelse(sex == "F", "1", "0"))

#select transformed variables only prior to imputation
#remove neutrophils (colinearity)
cov2 <- cov2 %>% dplyr::select(id, Sample_Sentrix_ID, age, sex, 
            Batch, logsmok, logbmi, logunits, simd2009v2rank, 
            Bcell, CD4T, CD8T, Eos, Mono, NK)
str(cov2)
cov2$sex <- as.factor(cov2$sex)

##6. Scale covariates: Scale numeric covariates
cov3 <- cov2 %>% modify_if(., is.numeric, scale)
sapply(cov3, function(y) sum(length(which(is.na(y)))))

#7. Impute missing data (KNN)
set.seed(1.234)
cov3 <- kNN(cov2, k = 5)
sapply(cov3, function(y) sum(length(which(is.na(y))))) #0


#remove imputation variables created
cov4 <- subset(cov3, select = id:NK)

#save

#create df for resid
cgdat <- cgdat %>% dplyr::filter(Sample_Sentrix_ID %in% cov4$Sample_Sentrix_ID)
dat <- dplyr::left_join(cgdat, cov4, by =  'id')
cov5 <- cov4 %>% dplyr::filter(Sample_Sentrix_ID %in% cgdat$Sample_Sentrix_ID)
sapply(dat, function(y) sum(length(which(is.na(y)))))

# #load kinship matrix
kin_model <- readRDS("kin_model_0901.rds")
kin_model2 <- kin_model * 2

#Residualise CpG
model <- lmekin(cg26561082 ~ age + sex + (1|id) + as.factor(Batch) + 
            logsmok + logbmi + logunits + simd2009v2rank + Bcell + 
            CD4T + CD8T + Eos + Mono + NK,
            varlist = kin_model2, 
            data = dat)

res <- residuals(model)
names(res) <- dat$Sample_Sentrix_ID.x
res <- data.frame(res)
res2 <- scale(res)
res2 <- data.frame(res2)
res2$Sample_Sentrix_ID <- row.names(res2)

#    # Protein Data
prot <- prot %>% dplyr::select(id, P19827)
prot$id <- as.character(prot$id)
# prot <- left_join(prot, ids, by = "id")
prot <- prot %>% dplyr::filter(id %in% dat$id)

#add cov 
protdat <- dplyr::left_join(prot, cov5, by = "id")
sapply(protdat, function(y) sum(length(which(is.na(y)))))


#Residualise protein
protmodel <- lmekin(P19827 ~ age + sex + (1|id) +
            logsmok + logbmi + logunits + simd2009v2rank,
            varlist = kin_model2, 
            data = protdat)
protres <- residuals(protmodel)
names(protres) <- protdat$Sample_Sentrix_ID
protres <- data.frame(protres)
protres2 <- scale(protres)
protres2 <- data.frame(protres2)
protres2$Sample_Sentrix_ID <- row.names(protres)

#Load pollutant data
poll <- fread("full_poll_2905.csv")
df <- dplyr::left_join(res2, poll, by = "Sample_Sentrix_ID")
df2 <- dplyr::left_join(protres2, df, by = "Sample_Sentrix_ID")
df2 <- df2 %>% dplyr::rename(prot_res = protres)

#2. Run models

med_res <- list()
pollutants <- c("PM25", "NO2", "NO", "PM10", "NO3_F", "O3")

# # 1. Mediator Model (Exposure -> Mediator)
for (pollutant in pollutants) {
formula_obj <- as.formula(paste0("res ~", pollutant))
model_mediator <- lm(formula_obj, data = df2)

# # 2. Outcome Model (Exposure + Mediator -> Outcome)
formula_obj <- as.formula(paste0("prot_res~", pollutant, " + res"))
model_outcome <- lm(formula_obj, data = df2)


set.seed(123) #for reproducibility of bootstrap results

model_mediator$call$formula <- formula(model_mediator)
model_outcome$call$formula <- formula(model_outcome)

mediation_results <- mediate(
  model.m = model_mediator,
  model.y = model_outcome,
  treat = paste0(pollutant),
  mediator = "res",
  boot = TRUE,
  sims = 5000 
)
med_res[[pollutant]] <- summary(mediation_results)
}

#outcome = Protein
#mediator = CpG
#Treatment = Pollutant

#save
saveRDS(med_res, file = "mediation_results.RDS")

#tabulate mediation results
tidy_mediation_table <- function(med_model) {
  # Bind the summary estimates into a tibble
  tibble(
    Effect = c("ACME (Average)", "ADE (Average)", "Total Effect", "Prop. Mediated"),
    Estimate = c(med_model$d.avg, med_model$z.avg, med_model$tau.coef, med_model$n.avg),
    Lower_CI = c(med_model$d.avg.ci[1], med_model$z.avg.ci[1], med_model$tau.ci[1], med_model$n.avg.ci[1]),
    Upper_CI = c(med_model$d.avg.ci[2], med_model$z.avg.ci[2], med_model$tau.ci[2], med_model$n.avg.ci[2]),
    p_value = c(med_model$d.avg.p, med_model$z.avg.p, med_model$tau.p, med_model$n.avg.p)
  )
}

res_df <- purrr::map(med_res, tidy_mediation_table)
res_df <- dplyr::bind_rows(res_df, .id = "pollutant")
fwrite(res_df, file = "mediation_results.csv")

res <- fread("mediation_results.csv")
prop <- res %>% dplyr::filter(Effect == "Prop. Mediated")


#Diagram
med_res <- readRDS("mediation_results.RDS")

# Force objects back to raw format
for (pollutant in names(med_res)) {
  class(med_res[[pollutant]]) <- "mediate"
}

pdf("all_mediation_plots.pdf", width = 8, height = 6)

for (pollutant in names(med_res)) {
  mediation:::plot.mediate(med_res[[pollutant]], 
       main = paste("Mediation Analysis:", pollutant))
}

# Close the device to finalize and save the file
dev.off()
