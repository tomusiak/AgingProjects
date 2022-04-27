#Example corrected MDP-Seq analysis of a dataset investigating leukemia.

#OPTIONAL: List of mitochondrial genes to correct for.
# Use this to exclude tRNA genes if desired.

keep_mitogenes <- c("RNR1","RNR2", "ND1","ND2","CO1","CO2","ATP8",
                    "ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")

#Loading in libraries.
library("BiocManager")
library("S4Vectors")
library("Biostrings")
library("SummarizedExperiment")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("Biobase")
library("AnnotationDbi")
library("GenomicFeatures")
library(dplyr)
library(ggplot2)
library("RColorBrewer")
library(Rsubread)
library(ggrepel)
library(circlize)

#Loading in auxiliary code. The primary code implementing the algorithm for MDP-Seq
# lives here.
setwd("/home/atom/Desktop/AgingProjects/Cohen/")
source("mdpseq_background_correction.R")

#Loading in auxiliary code. This has a number of "helper functions" that are useful
# across many pipelines & analyses.
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

#Relocating to directory where .fastq files and metadata live.
setwd("/home/atom/Desktop/Data/burridge2016")
#bam_list <- makeBAMS(".",FALSE) #Converts fastq to BAM.
#mitogene_counts <- getCountsMitochondrial(bam_list,FALSE) #Creates mitochondrial gene count matrix.
#mdp_counts <- getCountsMDP(bam_list,FALSE) #Creates MDP count matrix.
#write.csv(mitogene_counts,"mitogene_counts.csv")
#write.csv(mdp_counts,"mdp_counts.csv")


#Loads in count matrices and metadata. Appropriately converts them to factors and aligns
# column names.
mdp_counts <- data.frame(read.csv("mdp_counts.csv", row.names=1))
mitogene_counts <- data.frame(read.csv("mitogene_counts.csv",row.names=1))
burridge_samplesheet <- read.csv("burridge2016_samplesheet.csv")
burridge_samplesheet$doxo <- as.factor(burridge_samplesheet$doxo)
colnames(mdp_counts) <- burridge_samplesheet$sample
colnames(mitogene_counts) <- burridge_samplesheet$sample

#Performs a standard differential expression analysis for the mitochondrial gene.
# Not necessary, but can give a rough sense of changes in mitochondrial gene expression.
# Helps for sanity-checking background correction.
mitogene_dds <- DESeqDataSetFromMatrix(mitogene_counts,
                                       colData = burridge_samplesheet,
                                       design = ~ vulnerable+vulnerable:doxo)
mitogene_dds <- DESeq(mitogene_dds)
mitogene_res <- data.frame(results(mitogene_dds))
mitogene_vsd <- varianceStabilizingTransformation(mitogene_dds) ###mdp-seq script
mitogene_resOrdered <- mitogene_res[order(abs(mitogene_res$padj)),]

#Creates a volcano plot for mitochondrial gene expression changes.
mitogene_res$color=factor(case_when(mitogene_res$padj < .05 & abs(mitogene_res$log2FoldChange) >= .60 ~ "purple",
                                    (mitogene_res$padj < .05 & abs(mitogene_res$log2FoldChange) < .60) ~ "red",
                                    (mitogene_res$padj >= .05 & abs(mitogene_res$log2FoldChange) >= .60) ~ "blue",
                                    (mitogene_res$padj >= .05 & abs(mitogene_res$log2FoldChange) < .60) ~ "gray"))
mitogene_res$delabel <- NA
mitogene_res$delabel[mitogene_res$color=="purple"] <- rownames(mitogene_res)[mitogene_res$color=="purple"]
ggplot(data=mitogene_res, aes(x=log2FoldChange, y=-log10(padj), color=color, label=delabel)) + 
  geom_point() +
  xlim(-4,4) + ylim(0,22) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_text(nudge_x=.20) +
  scale_colour_identity() +
  geom_vline(xintercept = .6,linetype="dotted") +
  geom_vline(xintercept = -.6,linetype="dotted") +
  labs(title="MitoGene Volcano Plot for Leukemia")

#Reads in GTF files. These are being pulled from my local machine, but for convenience I have
# added them to the github drive as well, in "Useful Scripts."
mitogene_gtf <-
  read.delim(
    "/home/atom/Desktop/Data/mitochondria_db.gtf",
    header = FALSE,
    comment.char = "#"
  ) #Note: This mitochondrial gtf contains only mitochondrial genes.
mdp_gtf <-
  read.delim("/home/atom/Desktop/Data/mdp_atg_noATC.gtf", header = FALSE)
# A GTF of MDPs. Obtained from Brendan. 
# Likely best to remove the mitochondrial gene GTFs from this.

#This next segment of code performs the background correction.
# generateEncompassTable() is not dependent on the data set. It is responsible for creating a 
# table that maps out which MDPs are located on which mitochondrial gene (hereafter, "mitogene")
# and by roughly what proportion.
encompass_table <- generateEncompassTable(mdp_gtf,mitogene_gtf)

#OPTIONAL: Remove some of the mitochondrial genes. Can be used for excluding tRNAs.
encompass_table <- data.frame(encompass_table[encompass_table$mitogene %in% keep_mitogenes,])
mitogene_counts <- data.frame(mitogene_counts[rownames(mitogene_counts) %in% keep_mitogenes,])

#Creates a table using the mitogene count matrix that maps out relative changes in mitogene
# expression across samples relative to the mean.
background_table <- determineBackgroundSignal(mitogene_counts)

#Performs background correction by removing the background signal (calculated above) from
# the MDP counts relative to their overlap in the mitochondrial gene (calculate in 
# generateEncompassTable()).
corrected_counts <- data.frame(performBackgroundCorrection(background_table,mdp_counts,encompass_table))

#Performs a differential gene expression analysis using the corrected counts.
mdp_dds_corrected <- DESeqDataSetFromMatrix(corrected_counts,
                                            colData = burridge_samplesheet,
                                            design = ~patient + doxo)
mdp_dds_corrected <- DESeq(mdp_dds_corrected)
mdp_res_corrected <-data.frame(results(mdp_dds_corrected))
mdp_res_corrected <- mdp_res_corrected[!is.na(mdp_res_corrected$padj),]
mdp_res_corrected$condition <- "Doxo Vulnerability"
write.csv(mdp_res_corrected,"DoxoburridgeMDPSeq.csv")
mdp_vsd_corrected <- varianceStabilizingTransformation(mdp_dds_corrected) ###mdp-seq script
mdp_resOrdered_corrected <- mdp_res_corrected[order(abs(mdp_res_corrected$log2FoldChange),decreasing=TRUE),]
mdp_top_corrected <- head(mdp_resOrdered_corrected, 10)

#Creates a volcano plot using the corrected counts.
mdp_res_corrected$color=factor(case_when(mdp_res_corrected$padj < .05 & abs(mdp_res_corrected$log2FoldChange) >= .60 ~ "purple",
                                         (mdp_res_corrected$padj < .05 & abs(mdp_res_corrected$log2FoldChange) < .60) ~ "red",
                                         (mdp_res_corrected$padj >= .05 & abs(mdp_res_corrected$log2FoldChange) >= .60) ~ "blue",
                                         (mdp_res_corrected$padj >= .05 & abs(mdp_res_corrected$log2FoldChange) < .60) ~ "gray"))
mdp_res_corrected$delabel <- NA
mdp_res_corrected$delabel[mdp_res_corrected$color=="purple"] <- rownames(mdp_res_corrected)[mdp_res_corrected$color=="purple"]
ggplot(data=mdp_res_corrected, aes(x=log2FoldChange, y=-log10(padj), color=color, label=delabel)) + 
  geom_point() +
  xlim(-5,5) + ylim(0,22) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .6,linetype="dotted") +
  geom_vline(xintercept = -.6,linetype="dotted") +
  geom_text(nudge_x=.2) +
  scale_colour_identity() +
  labs(title="MDPSeq Volcano Plot for Effect of Doxorubicin")

#This code performs an uncorrected MDPSeq.
mdp_dds_uncorrected <- DESeqDataSetFromMatrix(mdp_counts,
                                              colData = burridge_samplesheet,
                                              design = ~patient+doxo)
mdp_dds_uncorrected <- DESeq(mdp_dds_uncorrected)
mdp_res_uncorrected <-data.frame(results(mdp_dds_uncorrected))
mdp_res_uncorrected <- mdp_res_uncorrected[!is.na(mdp_res_uncorrected$padj),]
mdp_vsd_uncorrected <- varianceStabilizingTransformation(mdp_dds_uncorrected) ###mdp-seq script
mdp_resOrdered_uncorrected <- mdp_res_uncorrected[order(abs(mdp_res_uncorrected$padj)),]
mdp_top_uncorrected <- head(mdp_resOrdered_uncorrected, 10)
mdp_res_uncorrected$color=factor(case_when(mdp_res_uncorrected$padj < .05 & abs(mdp_res_uncorrected$log2FoldChange) >= .60 ~ "purple",
                                           (mdp_res_uncorrected$padj < .05 & abs(mdp_res_uncorrected$log2FoldChange) < .60) ~ "red",
                                           (mdp_res_uncorrected$padj > .05 & abs(mdp_res_uncorrected$log2FoldChange) >= .60) ~ "blue",
                                           (mdp_res_uncorrected$padj > .05 & abs(mdp_res_uncorrected$log2FoldChange) < .60) ~ "gray"))
mdp_res_uncorrected$delabel <- NA
mdp_res_uncorrected$delabel[mdp_res_uncorrected$color=="purple"] <- rownames(mdp_res_uncorrected)[mdp_res_uncorrected$color=="purple"]
ggplot(data=mdp_res_uncorrected, aes(x=log2FoldChange, y=-log10(padj), color=color, label=delabel)) + 
  geom_point() +
  xlim(-5,5) + ylim(0,22) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .6,linetype="dotted") +
  geom_vline(xintercept = -.6,linetype="dotted") +
  geom_text(nudge_x=.2) +
  scale_colour_identity() +
  labs(title="MDPSeq Volcano Plot for Leukemia - Uncorrected")


