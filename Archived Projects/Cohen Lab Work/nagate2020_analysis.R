#keep_mitogenes <- c("RNR1","RNR2", "ND1","ND2","CO1","CO2","ATP8",
#                    "ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")

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

#Loading in auxiliary code.
setwd("/home/atom/Desktop/AgingProjects/Cohen/")
source("mdpseq_background_correction.R")

setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")
setwd("/home/atom/Desktop/Data/nagate2020")
#bam_list <- makeBAMS(".",FALSE)
#mitogene_counts <- getCountsMitochondrial(bam_list,FALSE)
#mdp_counts <- getCountsMDP(bam_list,FALSE)
#write.csv(mitogene_counts,"mitogene_counts.csv")
#write.csv(mdp_counts,"mdp_counts.csv")

#deleteFASTQs()
#deleteBAMs()

mdp_counts <- read.csv("mdp_counts.csv", row.names=1)
mitogene_counts <- read.csv("mitogene_counts.csv",row.names=1)
nagate_samplesheet <- read.csv("nagate_samplesheet.csv")
colnames(mdp_counts) <- nagate_samplesheet$sample
colnames(mitogene_counts) <- nagate_samplesheet$sample
nagate_samplesheet$leukemia <- as.factor(nagate_samplesheet$leukemia)
nagate_samplesheet$patient<- as.factor(nagate_samplesheet$patient)

mitogene_dds <- DESeqDataSetFromMatrix(mitogene_counts,
                                       colData = nagate_samplesheet,
                                       design = ~ patient + leukemia)
mitogene_dds <- DESeq(mitogene_dds)
mitogene_res <-data.frame(results(mitogene_dds))
mitogene_vsd <- varianceStabilizingTransformation(mitogene_dds) ###mdp-seq script
mitogene_resOrdered <- mitogene_res[order(abs(mitogene_res$padj)),]

mitogene_res$color=factor(case_when(mitogene_res$padj < .05 & abs(mitogene_res$log2FoldChange) >= .60 ~ "purple",
                                    (mitogene_res$padj < .05 & abs(mitogene_res$log2FoldChange) < .60) ~ "red",
                                    (mitogene_res$padj >= .05 & abs(mitogene_res$log2FoldChange) >= .60) ~ "blue",
                                    (mitogene_res$padj >= .05 & abs(mitogene_res$log2FoldChange) < .60) ~ "gray"))
mitogene_res$delabel <- NA
mitogene_res$delabel[mitogene_res$color=="purple"] <- rownames(mitogene_res)[mitogene_res$color=="purple"]
ggplot(data=mitogene_res, aes(x=log2FoldChange, y=-log10(padj), color=color, label=delabel)) + 
  geom_point() +
  xlim(-3,3) + ylim(0,22) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_text(nudge_x=.15) +
  scale_colour_identity() +
  geom_vline(xintercept = .6,linetype="dotted") +
  geom_vline(xintercept = -.6,linetype="dotted") +
  labs(title="MitoGene Volcano Plot for Lung Cancer- Resident T Cells")

mitogene_gtf <-
  read.delim(
    "/home/atom/Desktop/Data/mitochondria_db.gtf",
    header = FALSE,
    comment.char = "#"
  )
mdp_gtf <-
  read.delim("/home/atom/Desktop/Data/mdp_atg_noATC.gtf", header = FALSE)

encompass_table <- generateEncompassTable(mdp_gtf,mitogene_gtf)
#encompass_table <- encompass_table[encompass_table$mitogene %in% keep_mitogenes,]
background_table <- determineBackgroundSignal(mitogene_counts)
corrected_counts <- data.frame(performBackgroundCorrection(background_table,mdp_counts,encompass_table))

mdp_dds_corrected <- DESeqDataSetFromMatrix(corrected_counts,
                                            colData = nagate_samplesheet,
                                            design = ~ patient + leukemia)
mdp_dds_corrected <- DESeq(mdp_dds_corrected)
mdp_res_corrected <-data.frame(results(mdp_dds_corrected))
mdp_res_corrected <- mdp_res_corrected[!is.na(mdp_res_corrected$padj),]
mdp_res_corrected$condition <- "Leukemia"
write.csv(mdp_res_corrected,"LeukemiaNagateMDPSeq.csv")
mdp_vsd_corrected <- varianceStabilizingTransformation(mdp_dds_corrected) ###mdp-seq script
mdp_resOrdered_corrected <- mdp_res_corrected[order(abs(mdp_res_corrected$padj)),]
mdp_top_corrected <- head(mdp_resOrdered_corrected, 10)
mdp_res_corrected$color=factor(case_when(mdp_res_corrected$padj < .05 & abs(mdp_res_corrected$log2FoldChange) >= .60 ~ "purple",
                                         (mdp_res_corrected$padj < .05 & abs(mdp_res_corrected$log2FoldChange) < .60) ~ "red",
                                         (mdp_res_corrected$padj >= .05 & abs(mdp_res_corrected$log2FoldChange) >= .60) ~ "blue",
                                         (mdp_res_corrected$padj >= .05 & abs(mdp_res_corrected$log2FoldChange) < .60) ~ "gray"))
mdp_res_corrected$delabel <- NA
mdp_res_corrected$delabel[mdp_res_corrected$color=="purple"] <- rownames(mdp_res_corrected)[mdp_res_corrected$color=="purple"]
ggplot(data=mdp_res_corrected, aes(x=log2FoldChange, y=-log10(padj), color=color, label=delabel)) + 
  geom_point() +
  xlim(-4,4) + ylim(0,22) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .6,linetype="dotted") +
  geom_vline(xintercept = -.6,linetype="dotted") +
  geom_text(nudge_x=.2) +
  scale_colour_identity() +
  labs(title="MDPSeq Volcano Plot for (Nagate) Leukemia PBMCs")

mdp_dds_uncorrected <- DESeqDataSetFromMatrix(mdp_counts,
                                              colData = nagate_samplesheet,
                                              design = ~ patient + leukemia)
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
  xlim(-4,4) + ylim(0,22) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .6,linetype="dotted") +
  geom_vline(xintercept = -.6,linetype="dotted") +
  geom_text(nudge_x=.2) +
  scale_colour_identity() +
  labs(title="MDPSeq Volcano Plot for (Nagate) Leukemia PBMCs- Uncorrected")

#mess w/ mass spec data
ms_database <- read.csv("~/Desktop/Data/ms_database.csv", comment.char="#")
top_results <- mdp_res_corrected[mdp_res_corrected$padj < .1 & abs(mdp_res_corrected$log2FoldChange>.6),]
top_results$peptide <- rownames(top_results)
ms_database$mdp_id <- gsub(">Peptide", "", ms_database$peptide)
joined_database <- merge(ms_database,top_results)

top_hits_nagate <- top_results
