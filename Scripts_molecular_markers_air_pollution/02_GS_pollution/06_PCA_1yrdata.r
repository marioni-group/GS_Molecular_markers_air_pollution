#PCA of pollutants


#use micromamba environment: in terminal run:

eval "$(micromamba shell hook --shell bash)"
micromamba activate R
R

# if (!requireNamespace('BiocManager', quietly = TRUE))
# install.packages('BiocManager')
# BiocManager::install('PCAtools')

library(tidyverse)
library(data.table)
library(PCAtools)

#Load data
pollution <- fread("annual_avg_exp_foranalysis_scaled_0425.csv")

#set up data
# pollution <- pollution %>% dplyr::select(-year_range)
df1 <- data.table::transpose(pollution, keep.names = "pollutant", make.names = "ID")
df1 <- as.data.frame(df1)

row.names(df1) <- df1$pollutant
poll_list <- row.names(df1)

df2 <- df1[,2:22072]
pmetadata = data.frame("M" = rep(1, dim(df2)[2]), row.names = colnames(df2))

#carry out PCA analysis
pca_prot <- pca(df2, scale = TRUE, center = TRUE, metadata = pmetadata)
pca_summary <- summary(pca_prot)
prop_var <- round(pca_prot$variance, digits = 4)
prop_var
head(prop_var)
head(pca_prot$loadings)

pca_loadings <- pca_prot$loadings
pca_loadings$pollutant <- row.names(pca_loadings)

#save out
saveRDS(pca_prot, file = "pollution_pca.RDS")
fwrite(pca_loadings, file = "pca_loadings.csv")

prop_var_df <- data.frame(
    PC = names(prop_var),
    proportion = as.numeric(prop_var)
)
fwrite(prop_var_df, file ="pca_variance_prop.csv")


pdf("poll_screeplot.pdf")
screeplot(pca_prot, 
          drawCumulativeSumLine = FALSE, drawCumulativeSumPoints = FALSE) +
    geom_line(aes(x = 1:length(pca_prot$components), y = as.numeric(pca_prot$variance)))
dev.off()

pdf("poll_biplot.pdf")
biplot(pca_prot,
                showLoadings = TRUE, 
                ntopLoadings = 10,
                lab = NULL)
dev.off()

percent_load <- getComponents(pca_prot)
