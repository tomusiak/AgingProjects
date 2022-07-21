#This code analyzes data obtained from the Illumina EPIC array chip on four different CD8+ T cell
# subsets - Naive, CM, EM, and TEMRA cells. It identifies CpGs that are associated with T cell differentiation
# and removes them. It then filters on remaining CpGs that are associated with aging. In so doing, the program
# aims to create a list of CpGs that can be used for an epigenetic aging clock that are not associated with T
# cell differentiation.

#Grabs some useful scripts.
source("AgingProjects/Useful Scripts/generally_useful.R")

#Sets location of data 
setwd("Data/") #Sets directory.
#Libraries to import.
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(plyr)
library(reshape2)
library(WGCNA)
library(limma)
library(dplyr)
library("RColorBrewer")
library(readxl)
library(ggplot2)
library(readr)
library(umap)
require(clusterExperiment)
library(stringr)
library(tidyr)
library(Rfast)
library(minfi)

#Reading in database data for training and testing.
complete_cpg_table <- data.table::fread("ClockConstruction/complete_cpg_table.csv",header=TRUE) 
row.names(complete_cpg_table) <- complete_cpg_table$cpg
filtered_cpg_table <- complete_cpg_table
filtered_cpg_table <- filtered_cpg_table[,-c(1,2)]
complete_sample_table <- read.csv("ClockConstruction/complete_sample_table.csv",row.names=1)
filtered_cpg_table[filtered_cpg_table == "NULL"] <- NA
filtered_cpg_table[filtered_cpg_table == "null"] <- NA

#Read in clock data and metadata. Merge them and process them.
clock_data <- read.csv("subset_data/clock_data.csv")
subset_metadata <- read.csv("subset_data/subset_metadata.csv")
all_data <- merge(clock_data,subset_metadata)
all_data$type <- as.character(all_data$type)
all_data$type <- factor(all_data$type, levels=c("naive", "central_memory", "effector_memory","temra"))
beta_values <- read.csv("subset_data/beta_values.csv", row.names=1)
#m_values <- read.csv("subset_data/m_values.csv", row.names=1)

#Filter out unwanted data from metadata, mapping, and beta values.
keep <- all_data$SampleID[all_data$SampleID != "D3" & all_data$sabgal_sample == FALSE]

#Note - D3 was a mistakenly pipetted sample. SAbGal-high samples were originally included to
# determine if changes were occurring on a DNA methylation level, but ultimately only one
# SA-bGal sample was included and thus rigorous statistics cannot be performed.
all_data <- all_data[all_data$SampleID %in% keep,]
beta_values <- beta_values[,colnames(beta_values) %in% keep]
beta_values <- beta_values[,order(all_data$type)]
#m_values <- m_values[,colnames(m_values) %in% keep]
#m_values <- m_values[order(all_data$type)]
all_data <- all_data[order(all_data$type),]

#Here we will do differential methylation analysis on subsets, but subset down to CpGs that
# in the 450k panel.
validation_sample_table <- read.csv("ClockConstruction/validation_sample_table.csv",row.names=1)
validation_cpg_table <- read.csv("ClockConstruction/validation_cpg_table.csv",row.names=1)
validation_cpgs <- rownames(validation_cpg_table)
beta_values <- beta_values[rownames(beta_values) %in% rownames(validation_cpg_table),]

#Going to perform UMAP analysis on overall data set to determine where most variation is
# found. Intuitively, we will identify it in the T cell subsets that we assessed.
beta_rotated <- t(beta_values)
umap <- umap(beta_rotated,random_state=42)
umap_plot_df <- data.frame(umap$layout) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::inner_join(all_data, by = "SampleID")
umap_plot_df$type <- as.character(umap_plot_df$type)
umap_plot_df$type <- factor(umap_plot_df$type, levels=c("naive", "central_memory", 
                                                        "effector_memory","temra"))
ggplot(
  umap_plot_df,
  aes(x = X1, y = X2, color=type)) +
  labs(x="UMAP Component 1", y="UMAP Component 2", title = "UMAP Visualization - Before Filtering") +
  theme_classic() +
  geom_point(size=8) # Plot individual points to make a scatterplot
ggplot(umap_plot_df,aes(x=type,y=X1, color=type)) + 
  theme_classic() +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point(size=8)

#In the figures above, we found that UMAP component 2 seems to correlate with cell
# differentiation. We can use this to create a "pseudotime analysis" trajectory of bulk
# DNA methylation data where we look at how CpGs are differentially methylated with
# differentiation. We can then perform differential methylation analysis on the resulting
# sites to assess changes.
all_data$differentiation <- umap_plot_df$X1
celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45804","R45504",
                                              "R45553","R45741","R45805"))
diff_group <- all_data$differentiation
design <- model.matrix(~0 + donor_group + diff_group)
fit.reduced <- lmFit(beta_values,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
diff_exp <-topTable(fit.reduced,coef=8,number=1000000)

#Now we can remove all CpGs that change with differentiation.
changed_cpgs <- diff_exp[diff_exp$adj.P.Val < .1, ]
`%!in%` <- Negate(`%in%`)
filtered_cpgs <- beta_values[(rownames(beta_values) %!in% rownames(changed_cpgs)),]

#What if we now find best CpGs that change with age?
# Note - we can remove this step, and it may be best to as using this method inherently biases for 
# CpGs that specifically change with immune aging (although not with immune differentiation!)
celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45804","R45504",
                                              "R45553","R45741","R45805"))
age_design <- model.matrix(~0 + celltype_group + all_data$age)
fit_aging <- lmFit(filtered_cpgs,age_design)
fit_aging <- eBayes(fit_aging, robust=TRUE)
summary(decideTests(fit_aging))
diff_exp_aging <-topTable(fit_aging,coef=5,number=100000)
diff_exp_aging_order <- diff_exp_aging[order(diff_exp_aging$adj.P.Val),]
age_cpgs <- diff_exp_aging_order[diff_exp_aging_order$adj.P.Val < .2, ]
filtered_cpgs <- filtered_cpgs[rownames(filtered_cpgs) %in% rownames(age_cpgs),]

#Quickly double-checking if UMAP now segregates cell types..
nodiff_beta_rotated <- t(filtered_cpgs)
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
  labs(x="UMAP Component 1", y="UMAP Component 2", title = "UMAP Visualization - After Filtering") +
  theme_classic() +
  geom_point(size=8) # Plot individual points to make a scatterplot
ggplot(nodiff_umap_plot_df,aes(x=X1,y=X2, color=type)) + 
  theme_classic() +
  labs(x="UMAP Component 2", y="Age",title="UMAP Component 2 Tracks Age") +
  geom_point()
ggplot(nodiff_umap_plot_df,aes(x=age,y=X1)) + 
  theme_classic() +
  labs(x="Age", y="UMAP Component 1",title="UMAP Component 1 Tracks Age") +
  geom_point()

# ml_cpg_table <- ml_cpg_table[rownames(ml_cpg_table) %in% filtered_cpgs,]

#Writes list of candidate CpGs to be used for new clock construction.
write.csv(rownames(filtered_cpgs),"ClockConstruction/nodiff_cpgs.csv")

