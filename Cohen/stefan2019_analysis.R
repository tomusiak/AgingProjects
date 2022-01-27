keep_mitogenes <- c("RNR1","RNR2", "ND1","ND2","CO1","CO2","ATP8",
                    "ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")

library("BiocManager")
library("S4Vectors")
library("IRanges")
library("GenomeInfoDb")
library("GenomicRanges")
library("BiocGenerics")
library("BiocParallel")
library("XVector")
library("Biostrings")
library("Rsamtools")
library("matrixStats")
library("DelayedArray")
library("SummarizedExperiment")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("Biobase")
library("AnnotationDbi")
library("GenomicFeatures")
library(dplyr)
library("pheatmap")
library(ggplot2)
library("RColorBrewer")
library(Rsubread)
library(ggrepel)
library(umap)

#Loading in auxiliary code.
setwd("/home/atom/Desktop/AgingProjects/Cohen/")
source("mdpseq_background_correction.R")

setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

setwd("/home/atom/Desktop/Data/stefan2019")
sample_sheet <- read.csv("sample_sheet.csv")
sample_sheet$donor <- as.factor(sample_sheet$donor)
bam_list <- makeBAMS(".")
#mitogene_counts <- getCountsMitochondrial(bam_list,FALSE)
mdp_counts <- getCountsMDP(bam_list,FALSE)

mitogene_count_matrix <- read.csv("~/Desktop/Data/stefan2019/mitogene_count_matrix.csv", row.names=1)
mitogene_count_matrix <- mitogene_count_matrix[rownames(mitogene_count_matrix) %in% keep_mitogenes,]
mitogene_dds <- DESeqDataSetFromMatrix(mitogene_count_matrix,
                                       colData = sample_sheet,
                                       design = ~ donor + ifna)
mitogene_dds <- DESeq(mitogene_dds)
mitogene_res <-results(mitogene_dds)
mitogene_vsd <- varianceStabilizingTransformation(mitogene_dds) ###mdp-seq script
mitogene_resOrdered <- mitogene_res[order(abs(mitogene_res$padj)),]

cols <- densCols(mitogene_res$log2FoldChange, -log10(mitogene_res$padj),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= mitogene_res$log2FoldChange, 
     y = -log10(mitogene_res$padj), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=mitogene_res$pch, cex=0.4)
gn.selected <- abs(mitogene_res$log2FoldChange) >.5 & mitogene_res$padj < .05
text(mitogene_res$log2FoldChange[gn.selected],
     -log10(mitogene_res$padj)[gn.selected],
     lab=rownames(mitogene_res)[gn.selected ], cex=0.6)

mitogene_gtf <-
  read.delim(
    "/home/atom/Desktop/Data/mitochondria_db.gtf",
    header = FALSE,
    comment.char = "#"
  )
mdp_gtf <-
  read.delim("/home/atom/Desktop/Data/mdp_atg_noATC.gtf", header = FALSE)

encompass_table <- generateEncompassTable(mdp_gtf,mitogene_gtf)
encompass_table <- encompass_table[encompass_table$mitogene %in% keep_mitogenes,]
background_table <- determineBackgroundSignal(mitogene_count_matrix)
corrected_counts <- data.matrix(performBackgroundCorrection(background_table,mdp_counts,encompass_table))

mdp_dds_corrected <- DESeqDataSetFromMatrix(corrected_counts,
                                            colData = sample_sheet,
                                            design = ~ donor + ifna)
mdp_dds_corrected <- DESeq(mdp_dds_corrected)
mdp_res_corrected <-results(mdp_dds_corrected)
mdp_vsd_corrected <- varianceStabilizingTransformation(mdp_dds_corrected) ###mdp-seq script
mdp_resOrdered_corrected <- mdp_res_corrected[order(abs(mdp_res_corrected$padj)),]
mdp_top_corrected <- head(mdp_resOrdered_corrected, 10)

cols <- densCols(mdp_res_corrected$log2FoldChange, -log10(mdp_res_corrected$padj),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= mdp_res_corrected$log2FoldChange, 
     y = -log10(mdp_res_corrected$padj), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-3,3),
     ylim=c(0,20),
     pch=mdp_res_corrected$pch, cex=0.4)
gn.selected <- abs(mdp_res_corrected$log2FoldChange) >.75 & mdp_res_corrected$padj < .05
text(mdp_res_corrected$log2FoldChange[gn.selected],
     -log10(mdp_res_corrected$padj)[gn.selected],
     lab=rownames(mdp_res_corrected)[gn.selected ], cex=0.6)

mdp_dds_uncorrected <- DESeqDataSetFromMatrix(mdp_counts,
                                            colData = sample_sheet,
                                            design = ~ donor + ifna)
mdp_dds_uncorrected <- DESeq(mdp_dds_uncorrected)
mdp_res_uncorrected <-results(mdp_dds_uncorrected)
mdp_vsd_uncorrected <- varianceStabilizingTransformation(mdp_dds_uncorrected) ###mdp-seq script
mdp_resOrdered_uncorrected <- mdp_res_uncorrected[order(abs(mdp_res_uncorrected$padj)),]
mdp_top_uncorrected <- head(mdp_resOrdered_uncorrected, 10)
