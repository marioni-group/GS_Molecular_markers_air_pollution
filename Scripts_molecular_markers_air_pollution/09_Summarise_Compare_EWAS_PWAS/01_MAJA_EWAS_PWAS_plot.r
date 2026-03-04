#Combine EWAS and PWAS plots together

library(tidyverse)
library(data.table)
library(utils)
library(here)
library(ComplexHeatmap)

#EWAS plot
results_full2 <- fread("/maja_fullresults_1307_newthres.csv")

#Data-prep
results_plot <- results_full2 %>%
  dplyr::group_by(CpG) %>%
  dplyr::filter(any(association2 == "yes")) %>%
  dplyr::ungroup()
length(unique(results_plot$CpG)) 

results_plot <- results_plot %>%
  dplyr::mutate(association_value = ifelse(
    association2 == "yes", Betas, NA
  ))
head(results_plot)

#Tidy up for plotting
results_plot <- results_plot %>% 
    select(-c(association, lower, upper)) %>%
    mutate(Trait = ifelse(Trait == "logsmok", "Smoking", Trait))

#add Gene names for CpGs
anno <- fread("/majaewas_highpip_anno_0812_newthresh.csv")
head(anno)
anno <- anno %>% select(CpG, UCSC_RefGene_Name)
results_plot2 <- dplyr::left_join(results_plot, anno, by = "CpG", relationship = "many-to-many")
results_plot2 <- distinct(results_plot2)
results_plot2$UCSC_RefGene_Name <- str_remove_all(results_plot2$UCSC_RefGene_Name, ";.*$") #strip repeated gene names
results_plot2 <- results_plot2 %>% 
  mutate(UCSC_RefGene_Name = ifelse(UCSC_RefGene_Name == "", "nil", UCSC_RefGene_Name),
    plot_label = paste(CpG, UCSC_RefGene_Name, sep = ":"))


#Visual matrix with NA values
rmv <- results_plot2 %>%
  dplyr::select(plot_label, Trait, association_value) %>%
  tidyr::pivot_wider(
    names_from = Trait,
    values_from = association_value #where NA values are used in place of betas where no significant assocaition found
  )

#Create matrix of beta values including NAs 

rmv2 <- rmv %>% dplyr::select(-plot_label)
dims <- dim(rmv2)
rmv2 <- as.matrix(rmv2)
rmv2 <- as.numeric(rmv2)
dim(rmv2) <- dims 

#add protein names and pollutant names
rownames(rmv2) <- rmv$plot_label
colnames(rmv2) <- colnames(rmv)[2:10]

#Value matrix without NAs
rmi <- results_plot2 %>%
  dplyr::select(plot_label, Trait, Betas) %>%
  tidyr::pivot_wider(
    names_from = Trait,
    values_from = Betas
  )

#Add in protein names  
#Create matrix of beta values
rmi1 <- rmi %>% dplyr::select(-plot_label)
dims <- dim(rmi1)
rmi2 <- as.matrix(rmi1)
rmi2 <- as.numeric(rmi2)
dim(rmi2) <- dims 

#check orders of rmi and rmv are the same
identical(rmi$plot_label, rmv$plot_label) #TRUE

#add protein names and pollutant names
rownames(rmi2) <- rmi$plot_label
colnames(rmi2) <- colnames(rmi)[2:10]

#transpose
rmi2 <- t(rmi2)
rmv2 <- t(rmv2)
#Heatmap

#set row order
orderrows <- c("O3", "NO", "NO2", "SO2", "PM10", "NO3_F", "NO3_C", "PM25", "Smoking")
rmi2 <- rmi2[orderrows,]
rmv2 <- rmv2[orderrows,]

H1 <- ComplexHeatmap::Heatmap(rmi2, name = "Beta value",
               column_title = "CpGs",
               column_title_gp = gpar(fontsize = 18),
               row_title = "",
               col = circlize::colorRamp2(c(-0.2,0,0.2), c("blue", "white", "red")),
               row_names_gp = grid::gpar(fontsize = 18),
               column_names_gp = grid::gpar(fontsize = 12),
               show_column_dend = FALSE,
               row_names_side = "left",
               row_dend_side = "left",
               cluster_rows = FALSE,
              #  show_heatmap_legend = FALSE,
           cell_fun = function(j, i, x, y, width, height, fill) {
            if (is.na(rmv2[i,j])) {
            grid::grid.rect(x = x, y = y, width = width, height = height, 
                      gp = grid::gpar(fill = "grey90", col = NA))
          }
        })
 ComplexHeatmap::draw(H1)
ht1_drawn <- draw(H1)
clustered_order <- row_order(ht1_drawn)
clustered_rownames <- rownames(rmi2)[clustered_order]



#Load results
results <- fread("proteins_pollution/maja_fullresults_scpoll_burninadj_0812_newthresh.csv")
results <- results %>%
  mutate(associated_direction = ifelse(Betas < 0 & association == "yes", "Negative",
                                        ifelse(Betas > 0 & association == "yes", "Positive", NA)))

#Heatmaps
#Create two matrices - first for visualisation where NAs replace beta values for non-highpip associations
# the second for clustering

results_plot <- results %>%
  group_by(protein) %>%
  dplyr::filter(any(association2 == "yes")) %>%
  ungroup()
length(unique(results_plot$protein)) #45 proteins retained

results_plot <- results_plot %>%
  select(-c(association, lower, upper)) %>%
  mutate(association_value = ifelse(
    association2 == "yes", Betas, NA
  ))

head(results_plot)
results_plot <- results_plot %>% mutate(Trait = ifelse(Trait == "logsmok", "Smoking", Trait))
#Data frame with all results for any protein with a significant association:
results_plot #(see above)
#results_plot

#Add in protein names  
short_annots <- fread("groups_idmapping_01_25_labelled_simplified.csv") #
short_annots <- short_annots %>% dplyr::rename(protein = original_id)
short_annots <- short_annots %>% dplyr::filter(protein %in% results_plot$protein)
short_annots2 <- short_annots %>% dplyr::select(protein, gene_names_simpl)
short_annots2$gene_names_simpl <- str_remove_all(short_annots2$gene_names_simpl, ";*")

results_plot2 <- left_join(results_plot, short_annots2, by = "protein")
results_plot2 <- distinct(results_plot2)

#add in protein categories
library(readr)
cat <- read_csv("unique_proteins_function.csv")
cat <- cat %>% select("UniProt ID", "simplified_function") %>%
  dplyr::rename(protein = "UniProt ID",
  category = "simplified_function")
cat <- cat %>% filter(protein %in% results_plot$protein)
length(unique(cat$category)) 
unique(cat$category)

#fix error
cat <- cat %>% mutate(category = ifelse(category == "vascular homeostatis", "vascular homeostasis", category))

#add categories to results_plot
results_plot2 <- left_join(results_plot2, cat, by = "protein")
results_plot2$category <- factor(results_plot2$category, levels = c("immune response", "complement system", "coagulation", "metabolism", "endocrine", "gas transport", "vascular homeostasis", "extra-cellular matrix"))
results_plot2 <- results_plot2 %>% dplyr::rename(fntn = category)
results_plot2 <- results_plot2 %>% arrange(fntn)
results_plot2 <- results_plot2 %>% mutate(gene_names_simpl = ifelse(gene_names_simpl == "", "Ig-like", gene_names_simpl))

fwrite(results_plot2, file = "maja_results_categorisedprot_0812_newthresh.csv")

results_plot2 %>% 
   group_by(fntn) %>%
   summarise(n = length(unique(protein)))

#Visual matrix with NA values
rmv <- results_plot2 %>%
  dplyr::select(protein, Trait, association_value) %>%
  pivot_wider(
    names_from = Trait,
    values_from = association_value #where NA values are used in place of betas where no significant association found
  )

#Create matrix of beta values including NAs 
short_annots2 <- short_annots2 %>% mutate(gene_names_simpl = ifelse(gene_names_simpl == "", "Ig-like", gene_names_simpl))
rmv1 <- left_join(short_annots2, rmv, by = "protein")
rmv2 <- rmv1 %>% dplyr::select(-c(protein, gene_names_simpl))
dims <- dim(rmv2)
rmv2 <- as.matrix(rmv2)
rmv2 <- as.numeric(rmv2)
dim(rmv2) <- dims 

#add protein names and pollutant names
rownames(rmv2) <- rmv1$gene_names_simpl
colnames(rmv2) <- colnames(rmv1)[3:11]

#Value matrix without NAs
rmi <- results_plot2 %>%
  dplyr::select(protein, Trait, Betas) %>%
   pivot_wider(
    names_from = Trait,
    values_from = Betas
  )

#Add in protein names  
#Create matrix of beta values
rmi <- left_join(short_annots2, rmi, by = "protein")
rmi1 <- rmi %>% dplyr::select(-c(protein, gene_names_simpl))
dims <- dim(rmi1)
rmi2 <- as.matrix(rmi1)
rmi2 <- as.numeric(rmi2)
dim(rmi2) <- dims 

#check orders of rmi and rmv are the same
identical(rmi$protein, rmv1$protein) #TRUE

#add protein names and pollutant names
rownames(rmi2) <- rmi$gene_names_simpl
colnames(rmi2) <- colnames(rmi1)

rmi2 <- t(rmi2)
rmv2 <- t(rmv2)

#match row order to cpg heatmap
rmi2 <- rmi2[orderrows,]
rmv2 <- rmv2[orderrows,]
identical(rownames(rmi2), rownames(rmv2))
identical(colnames(rmi2), colnames(rmv2))

#Heatmap
cat_order <- results_plot2 %>% select(gene_names_simpl, fntn)
cat_vec <- cat_order$fntn[match(colnames(rmi2), cat_order$gene_names_simpl)]
dend2 = ComplexHeatmap::cluster_within_group(rmi2, cat_vec)
library(paletteer)
clist <- paletteer_d("ggthemes::Classic_Color_Blind")
clist <- c("immune response" = "#006BA4FF", 
             "complement system" = "#5F9ED1FF",
             "coagulation" = "#C85200FF",
             "metabolism" = "#FFBC79FF",
             "endocrine" = "#FF800EFF",
             "gas transport" = "#ABABABFF",
             "vascular homeostasis" = "#898989FF",
             "extra-cellular matrix" = "#CFCFCFFF")

H2 <- Heatmap(rmi2, name = "Beta value",
               column_title = "Proteins",
               column_title_gp = gpar(fontsize = 18),
               row_title = "",
               col = circlize::colorRamp2(c(-0.2,0,0.2), c("blue", "white", "red")),
               row_names_gp = gpar(fontsize = 18),
               row_names_side = "left",
               cluster_rows = FALSE,
              cluster_columns = dend2, 
              column_split = 8,
               top_annotation = HeatmapAnnotation(
                "Protein functional category" = cat_vec, 
                col = list("Protein functional category" = clist),  
                annotation_legend_param = list(title = "Protein functional category",
                direction = "horizontal", nrow = 2, title_POSITION = "lefttop")),
               heatmap_legend_param = list(
    direction = "horizontal",
    legend_width = unit(4, "cm"),   # length of the horizontal bar
    title_position = "lefttop"      # optional
  ),
               column_names_gp = gpar(fontsize = 12),
               show_column_dend = FALSE,
           cell_fun = function(j, i, x, y, width, height, fill) {
            if (is.na(rmv2[i,j])) {
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(fill = "grey90", col = NA))
          }
        })
 draw(H2)

H1 + H2

pdf("/ewas_pwas_together5_newthresh.pdf", width = 20, height = 6, family = "Helvetica")
 draw(H1 + H2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE)
dev.off()

png("/ewas_pwas_together5_newthresh.png", width = 1400, height = 600, family = "Helvetica")
 draw(H1 + H2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = TRUE)
dev.off()



#Create multi-plot with variance explained

#load ewas variances
all <- fread("/ewas_variances_0725.csv")

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

ewas <- all %>%
  ggplot(aes(x = Measure, y = Variance, fill = Measure)) +
  scale_fill_manual(values = clist) +
  geom_col() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  theme_bw() +
  labs(x = "Pollutant")

#load protein variances
all <- fread("/maja_prot_variance_explained.csv")

#Plot
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

prot_plot  <- all %>%
  ggplot(aes(x = Measure, y = Variance, fill = Measure)) +
  scale_fill_manual(values = clist) +
  geom_col() +
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  theme_bw() +
  labs(x = "Pollutant")

#variance plots together
library(patchwork)
ewas + prot_plot


img <- png::readPNG("/ewas_pwas_together5_newthresh.png")
htmp_grob <- rasterGrob(img, interpolate = TRUE)
htmp <- wrap_elements(full = htmp_grob)



# # Combine with patchwork
p1 <- (ewas + theme(plot.margin = unit(c(25,20,40,25), "pt")) + prot_plot + theme(plot.margin = unit(c(5,20,40,25), "pt")))
p2 <- p1 /htmp + plot_layout(heights = c(1, 1.5))


png("/multi_panel_newthresh.png", width = 1400, height = 1000, family = "Helvetica")
p2 + plot_annotation(tag_levels = 'A')
dev.off()


#Plot for presentation 

protdf <- fread("/maja_results_categorisedprot_0812_newthresh.csv")
sigdf <- protdf %>% filter(association2 == "yes")
sigdf$Trait <- factor(sigdf$Trait, levels = c("Smoking", "NO3_C", "NO3_F", "PM25", "PM10", "O3", "NO", "NO2", "SO2")) 

clist <- paletteer_d("ggthemes::Classic_Color_Blind")
clist <- c("immune response" = "#006BA4FF", 
             "complement system" = "#5F9ED1FF",
             "coagulation" = "#C85200FF",
             "metabolism" = "#FFBC79FF",
             "endocrine" = "#FF800EFF",
             "gas transport" = "#ABABABFF",
             "vascular homeostasis" = "#898989FF",
             "extra-cellular matrix" = "#CFCFCFFF")


sigdf %>%
  ggplot(aes(x = Trait, fill = fntn)) +
  scale_fill_manual(values = clist) +
  geom_bar() + 
  theme_bw()
ggsave(file = "/prot_results_presplot.pdf")
