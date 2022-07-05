#Grabs some useful scripts.
source("AgingProjects/Useful Scripts/generally_useful.R")

#Sets location of data 
setwd("Data/") #Sets directory.
#Libraries to import.
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(methylGSA)
library(plyr)
library(reshape2)
library(WGCNA)
library(limma)
library(dplyr)
require("biomaRt")
library("RColorBrewer")
library(readxl)
library(ggplot2)
library(GOfuncR)
library(readr)
library(plotrix)
library(umap)
library(splines)
library(gam)
require(clusterExperiment)
library(stringr)
library(tidyr)
library(gplots)
library(minfi)

#Read in clock data and metadata. Merge them and process them.
clock_data <- read.csv("subset_data/clock_data.csv")
subset_metadata <- read.csv("subset_data/subset_metadata.csv")
all_data <- merge(clock_data,subset_metadata)
all_data$type <- as.character(all_data$type)
all_data$type <- factor(all_data$type, levels=c("naive", "central_memory", "effector_memory","temra"))
beta_values <- read.csv("subset_data/beta_values.csv", row.names=1)
m_values <- read.csv("subset_data/m_values.csv", row.names=1)

#Here we will do differential methylation analysis on subsets, but subset down to CpGs that
# in the 450k panel.
validation_sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
validation_cpg_table <- read.csv("ClockConstruction/cpg_table.csv",row.names=1)
validation_cpgs <- rownames(validation_cpg_table)
filtered_beta_values <- beta_values[rownames(beta_values) %in% rownames(validation_cpg_table),]

#Now we can remove all CpGs that change with differentiation.
changed_cpgs <- diff_exp[diff_exp$adj.P.Val < .05, ]
`%!in%` <- Negate(`%in%`)
filtered_cpgs <- filtered_beta_values[(rownames(filtered_beta_values) %!in% rownames(changed_cpgs)),]

#Quickly double-checking if UMAP now segregates cell types..
nodiff_beta_values <- beta_values[rownames(beta_values) %in% rownames(filtered_cpgs),]
nodiff_fit.reduced <- lmFit(nodiff_beta_values,validation_design)
nodiff_fit.reduced <- eBayes(nodiff_fit.reduced, robust=TRUE)
nodiff_beta_rotated <- t(nodiff_beta_values)
nodiff_umap <- umap(nodiff_beta_rotated)
nodiff_umap_plot_df <- data.frame(nodiff_umap$layout) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::inner_join(all_data, by = "SampleID")
nodiff_umap_plot_df$type <- as.character(nodiff_umap_plot_df$type)
nodiff_umap_plot_df$type <- factor(nodiff_umap_plot_df$type, 
                                   levels=c("naive", "central_memory", "effector_memory","temra"))

#Plotting all of this now.
ggplot(
  nodiff_umap_plot_df,
  aes(x = X1, y = X2, color=type )) +
  labs(x="UMAP Component 1", y="UMAP Component 2", title = "UMAP Visualization") +
  theme_classic() +
  geom_point() # Plot individual points to make a scatterplot

write.csv(rownames(nodiff_beta_values),"ClockConstruction/nodiff_cpgs.csv")