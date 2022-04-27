#Grabs some useful scripts.
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

#Sets location of data 
setwd("/home/atom/Desktop/Data/subset_data") #Sets directory.
#Libraries to import.
library(plyr)
library(reshape2)
library(limma)
library(dplyr)
library("RColorBrewer")
library(readxl)
library(ggplot2)
library(readr)
library(plotrix)
library(sesame)
library(umap)
library(splines)
library(gam)
require(clusterExperiment)
library(stringr)
library(tidyr)
library(rstatix)
library(ggpubr)
library(gplots)
library(minfi)

numtome <- read.csv("~/Desktop/Data/numtome.csv", comment.char="#")
cpg_annotation <- read.delim("~/Desktop/Data/subset_data/EPIC.hg38.manifest.tsv")
cpg_annotation <- cpg_annotation[,1:5]
desired_columns <- c(1,2,13,14,15,17)
numtome <- numtome[,desired_columns]
colnames(numtome) <- c("MDP","NUMT","chr","start","stop","strand")
numtome <- numtome[!duplicated(numtome$NUMT),]
numtome <- numtome[!duplicated(numtome$start),]
numtome <- numtome[!duplicated(numtome$stop),]
numtome <- numtome[,c("chr","MDP","start","stop","strand","NUMT")]
colnames(cpg_annotation) <- c("chr","start","stop","strand","probe_name")
cpg_table <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(cpg_table) <- c("cpg","gene")
upstream <- 500
downstream <- 2000
for (cpg in 1:nrow(cpg_annotation)) {
  cpg_chr <- cpg_annotation[cpg,1]
  cpg_start <- cpg_annotation[cpg,2]
  if(!is.na(cpg_start)) {
    relevant_genes <- numtome[numtome$chr == cpg_chr,]
    associated_genes <- relevant_genes[(relevant_genes$start - cpg_start <= downstream &
                                          relevant_genes$start - cpg_start >= 0) |
                                         (cpg_start - relevant_genes$start <= upstream &
                                            cpg_start - relevant_genes$start >= 0),]
    if (nrow(associated_genes) != 0) {
      associated_genes$cpg <- cpg_annotation[cpg,5]
      added_table <- associated_genes[,c(6,7)]
      cpg_table <- rbind(cpg_table,added_table)
    }
  }
}

beta_values <- read.csv("~/Desktop/Data/subset_data/beta_values.csv", row.names=1)

#Filter out unwanted data from metadata, mapping, and beta values.

#Read in clock data and metadata. Merge them and process them.
clock_data <- read.csv("~/Desktop/Data/subset_data/clock_data.csv")
subset_metadata <- read.csv("~/Desktop/Data/subset_data/subset_metadata.csv")
all_data <- merge(clock_data,subset_metadata)
all_data$type <- as.character(all_data$type)
all_data$type <- factor(all_data$type, levels=c("naive", "central_memory", "effector_memory","temra"))

keep <- all_data$SampleID[all_data$SampleID != "D3" & all_data$sabgal_sample == FALSE]
all_data <- all_data[all_data$SampleID %in% keep,]
beta_values <- beta_values[,colnames(beta_values) %in% keep]
beta_values <- beta_values[order(all_data$type)]
all_data <- all_data[order(all_data$type),]

celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45804","R45504",
                                              "R45553","R45741","R45805"))
sabgal_sample <- factor(all_data$sabgal_sample,levels=c("TRUE","FALSE"))
old_group <- factor(all_data$old,levels=c(TRUE,FALSE))
diff_group <- all_data$differentiation
design <- model.matrix(~0 + donor_group + celltype_group)

#Let's look at aging!
all_data$older <- c("Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Older","Older","Older","Younger","Older")
age_group <- all_data$older
age_design <- model.matrix(~0 + celltype_group + age_group)

fit_aging <- lmFit(beta_values,age_design)
fit_aging <- eBayes(fit_aging, robust=TRUE)
diff_exp_aging <-topTable(fit_aging,coef=5,number=1000000)
numt_changes <- diff_exp_aging[rownames(diff_exp_aging) %in% cpg_table$cpg,]

combined_mapping_table <- cpg_table
combined_mapping_table$MDP <- numtome$MDP[match(cpg_table$NUMT,numtome$NUMT)]
numt_changes$MDP <- combined_mapping_table$MDP[match(rownames(numt_changes),
                                                     combined_mapping_table$cpg)]

numt_changes$color=factor(case_when(numt_changes$adj.P.Val < .05 & abs(numt_changes$logFC) >= .10 ~ "purple",
                                        (numt_changes$adj.P.Val < .05 & abs(numt_changes$logFC) < .10) ~ "red",
                                        (numt_changes$adj.P.Val >= .05 & abs(numt_changes$logFC) >= .10) ~ "blue",
                                        (numt_changes$adj.P.Val >= .05 & abs(numt_changes$logFC) < .10) ~ "gray"))
numt_changes$delabel <- NA
numt_changes$delabel[numt_changes$color=="purple"] <- numt_changes$MDP[numt_changes$color=="purple"]
ggplot(data=numt_changes, aes(x=logFC, y=-log10(adj.P.Val), color=color,label=delabel)) + 
  geom_point() +
  xlim(-.50,.50) + ylim(0,20) +
  theme_classic(base_size=14)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .1,linetype="dotted") +
  geom_vline(xintercept = -.1,linetype="dotted") +
  geom_text(y=2) +
  scale_colour_identity() +
  labs(title="Volcano Plot for NUMTs Associated with MDPs w/ Aging")

fit_diff <- lmFit(beta_values,design)
fit_diff <- eBayes(fit_diff, robust=TRUE)
diff_exp_diff <-topTable(fit_diff,coef=10,number=1000000)
numt_changes_diff <- diff_exp_diff[rownames(diff_exp_diff) %in% cpg_table$cpg,]

numt_changes_diff$MDP <- combined_mapping_table$MDP[match(rownames(numt_changes_diff),
                                                     combined_mapping_table$cpg)]
numt_changes_diff$NUMT <- combined_mapping_table$NUMT[match(rownames(numt_changes_diff),
                                                          combined_mapping_table$cpg)]
numt_changes_diff$MDP <- sub('.Peptide', '',numt_changes_diff$MDP)
numt_changes_diff$NUMT <- sub('numt_smorf_', '',numt_changes_diff$NUMT)
numt_changes_diff$color=factor(case_when(numt_changes_diff$adj.P.Val < .05 & abs(numt_changes_diff$logFC) >= .20 ~ "purple",
                                    (numt_changes_diff$adj.P.Val <= .05 & abs(numt_changes_diff$logFC) <= .20) ~ "red",
                                    (numt_changes_diff$adj.P.Val >= .05 & abs(numt_changes_diff$logFC) >= .20) ~ "blue",
                                    (numt_changes_diff$adj.P.Val >= .05 & abs(numt_changes_diff$logFC) < .20) ~ "gray"))
numt_changes_diff$delabel <- NA
numt_changes_diff$delabel[numt_changes_diff$color=="purple"] <- numt_changes_diff$MDP[numt_changes_diff$color=="purple"]
ggplot(data=numt_changes_diff, aes(x=logFC, y=-log10(adj.P.Val), color=color,label=delabel)) + 
  geom_point() +
  xlim(-.50,.50) + ylim(0,8) +
  theme_classic(base_size=11)  +
  geom_hline(yintercept = 1.28,linetype="dotted") +
  geom_vline(xintercept = .1,linetype="dotted") +
  geom_vline(xintercept = -.1,linetype="dotted") +
  geom_text(size=3.7,position=position_nudge(x=-.04)) +
  scale_colour_identity() +
  labs(title="Volcano Plot for NUMT CpGs Associated with MDPs w/ Immune Differentiation")

grabCPGs <- function(NUMT) {
  return(cpg_table[cpg_table$NUMT == NUMT,2])
}

all_NUMTs <- numtome$NUMT
NUMT <- all_NUMTs[1]
NUMT <- beta_values[grabCPGs(NUMT),]
NUMT <- data.frame(t(NUMT))
NUMT$mean <- rowMeans(NUMT)
complete_table <- data.frame(t(NUMT$mean))
colnames(complete_table) <- rownames(NUMT)
rownames(complete_table) <- all_NUMTs[1]

for (x in 2:length(all_NUMTs)) {
  NUMT <- all_NUMTs[x]
  NUMT <- beta_values[grabCPGs(NUMT),]
  NUMT <- data.frame(t(NUMT))
  NUMT$mean <- rowMeans(NUMT)
  NUMT_table <- data.frame(t(NUMT$mean))
  colnames(NUMT_table) <- rownames(NUMT)
  rownames(NUMT_table) <- all_NUMTs[x]
  complete_table <- rbind(complete_table,NUMT_table)
}

complete_table <- na.omit(complete_table)
complete_table$MDP <- combined_mapping_table$MDP[match(rownames(complete_table),
                                                     numtome$NUMT)]
complete_table <- na.omit(complete_table)
MDPs <- complete_table %>% group_by(MDP) %>% 
  summarise(across(everything(), mean))
MDPs <- data.frame(MDPs)
rownames(MDPs) <- MDPs[,1]
MDPs <- data.frame(MDPs)[,-1]

fit_aging <- lmFit(MDPs,age_design)
fit_aging <- eBayes(fit_aging, robust=TRUE)
diff_exp_aging <-topTable(fit_aging,coef=5,number=1000000)
numt_changes <- diff_exp_aging[rownames(diff_exp_aging) %in% cpg_table$cpg,]

fit_diff <- lmFit(MDPs,design)
fit_diff <- eBayes(fit_diff, robust=TRUE)
mdp_exp_diff <-topTable(fit_diff,coef=10,number=1000000)
mdp_exp_diff <- mdp_exp_diff[-7,]
mdp_exp_diff$color=factor(case_when(mdp_exp_diff$adj.P.Val < .05 & abs(mdp_exp_diff$logFC) >= .20 ~ "purple",
                                         (mdp_exp_diff$adj.P.Val <= .05 & abs(mdp_exp_diff$logFC) <= .20) ~ "red",
                                         (mdp_exp_diff$adj.P.Val >= .05 & abs(mdp_exp_diff$logFC) >= .20) ~ "blue",
                                         (mdp_exp_diff$adj.P.Val >= .05 & abs(mdp_exp_diff$logFC) < .20) ~ "gray"))
mdp_exp_diff$delabel <- NA
rownames(mdp_exp_diff) <- sub('.Peptide', '',rownames(mdp_exp_diff) )
mdp_exp_diff$delabel[mdp_exp_diff$color=="purple"] <- rownames(mdp_exp_diff)[mdp_exp_diff$color=="purple"]
ggplot(data=mdp_exp_diff, aes(x=logFC, y=-log10(adj.P.Val), color=color,label=delabel)) + 
  geom_point() +
  xlim(-.50,.50) + ylim(0,8) +
  theme_classic(base_size=11)  +
  geom_hline(yintercept = 1.28,linetype="dotted") +
  geom_vline(xintercept = .20,linetype="dotted") +
  geom_vline(xintercept = -.20,linetype="dotted") +
  geom_text(size=3.7,position=position_nudge(x=-.05,y=0.03)) +
  scale_colour_identity() +
  labs(title="Volcano Plot for MDPs (all NUMTs averaged) w/ Immune Differentiation")
