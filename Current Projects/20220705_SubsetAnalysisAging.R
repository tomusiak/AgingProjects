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
require("biomaRt")
library("RColorBrewer")
library(readxl)
library(ggplot2)
library(GOfuncR)
library(readr)
library(plotrix)
library(umap)
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

#Filter out unwanted data from metadata, mapping, and beta values.
keep <- all_data$SampleID[all_data$SampleID != "D3" & all_data$sabgal_sample == FALSE]

#Note - D3 was a mistakenly pipetted sample. SAbGal-high samples were originally included to
# determine if changes were occurring on a DNA methylation level, but ultimately only one
# SA-bGal sample was included and thus rigorous statistics cannot be performed.
all_data <- all_data[all_data$SampleID %in% keep,]
beta_values <- beta_values[,colnames(beta_values) %in% keep]
beta_values <- beta_values[,order(all_data$type)]
m_values <- m_values[,colnames(m_values) %in% keep]
m_values <- m_values[order(all_data$type)]
all_data <- all_data[order(all_data$type),]

#Creating a differential methylation analysis, but now focusing solely on age.
celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45804","R45504",
                                              "R45553","R45741","R45805"))
all_data$older <- c("Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Older","Older","Older","Younger","Older")
age_group <- all_data$older
age_design <- model.matrix(~0 + celltype_group + age_group)

#Investigating and importing a table of CpGs-associated with gene promoters.
complete_table <- read.csv("subset_data/complete_table.csv", row.names=1)
cpg_table <- read.csv("subset_data/cpg_table.csv", row.names=1)
fit.reduced_age <- lmFit(complete_table,age_design)
fit.reduced_age <- eBayes(fit.reduced_age, robust=TRUE)
summary(decideTests(fit.reduced_age))
diff_exp_age <-topTable(fit.reduced_age,coef=5,number=400000)
diff_exp_order_age <- diff_exp_age[order(diff_exp_age$adj.P.Val),]
diff_exp_order_age <-
  diff_exp_order_age[!grepl("ENSG",rownames(diff_exp_order_age)),]
top_age <- rownames(head(diff_exp_order_age,12))
top_list_age<-getDiffMethylationListAge(top_age,cpg_table)
top_list_age<-drop_na(melt(top_list_age))

#Plotting which genes change the most with age, investigating specifically DNA methylation via
# Violin plots.
ggplot(top_list_age,aes(x=older,y=value,fill=older)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Age Group") +
  ylab("Methylated %") +
  ggtitle("Most Differentially-Methylated Genes - Aging") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

#Creates a volcano plot using the corrected counts.
diff_exp_order_age$color=factor(case_when(diff_exp_order_age$adj.P.Val < .05 & abs(diff_exp_order_age$logFC) >= .1 ~ "purple",
                                          (diff_exp_order_age$adj.P.Val < .05 & abs(diff_exp_order_age$logFC) < .1) ~ "red",
                                          (diff_exp_order_age$adj.P.Val >= .05 & abs(diff_exp_order_age$logFC) >= .1) ~ "blue",
                                          (diff_exp_order_age$adj.P.Val >= .05 & abs(diff_exp_order_age$logFC) < .1) ~ "gray"))
diff_exp_order_age$delabel <- NA
diff_exp_order_age$delabel[diff_exp_order_age$color=="purple"] <- rownames(diff_exp_order_age)[diff_exp_order_age$color=="purple"]
ggplot(data=diff_exp_order_age, aes(x=logFC, y=-log10(adj.P.Val), color=color)) + 
  geom_point() +
  xlim(-.5,.5) + ylim(0,20) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .1,linetype="dotted") +
  geom_vline(xintercept = -.1,linetype="dotted") +
  scale_colour_identity() +
  labs(title="Volcano Plot for Aging")

#Performing differential methylation analysis, but this time focusing on CpGs rather than promoters of
# genes.
fit_aging <- lmFit(beta_values,age_design)
fit_aging <- eBayes(fit_aging, robust=TRUE)
summary(decideTests(fit_aging))
diff_exp_aging <-topTable(fit_aging,coef=5,number=100000)
diff_exp_aging_order <- diff_exp_aging[order(diff_exp_aging$adj.P.Val),]
pvals_aging <- as.vector(diff_exp_aging_order[,5])
names(pvals_aging) <- rownames(diff_exp_aging_order)

#Investigating changes specifically in naive T cells only.
naive_keep <- all_data$SampleID[all_data$type == "naive"]
naive_data <- all_data[all_data$SampleID %in% naive_keep,]
naive_beta_values <- beta_values[,colnames(beta_values) %in% naive_keep]
naive_beta_values <- naive_beta_values[order(naive_data$type)]
age_naive <- factor(naive_data$older)
donor_naive <- factor(naive_data$donor)
design <- model.matrix(~0 + age_naive)
naive_fit_reduced <- lmFit(naive_beta_values,design)
naive_fit_reduced <- eBayes(naive_fit_reduced, robust=TRUE)
summary(decideTests(naive_fit_reduced))
diff_exp_naive <-topTable(naive_fit_reduced,coef=1,number=100000)
