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
library(circlize)

#Loading in auxiliary code.
setwd("/home/atom/Desktop/AgingProjects/Cohen/")
source("mdpseq_background_correction.R")

setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

setwd("/home/atom/Desktop/Data/sun2022")
#bam_list <- makeBAMS(".",FALSE)
#mitogene_counts <- getCountsMitochondrial(bam_list,FALSE)
#mdp_counts <- getCountsMDP(bam_list,FALSE)
#write.csv(mitogene_counts,"mitogene_counts.csv")
#write.csv(mdp_counts,"mdp_counts.csv")

mdp_counts <- read.csv("mdp_counts.csv", row.names=1)
mitogene_counts <- read.csv("mitogene_counts.csv",row.names=1)
sun_samplesheet <- read.csv("sun_samplesheet.csv")
colnames(mdp_counts) <- sun_samplesheet$sample
colnames(mitogene_counts) <- sun_samplesheet$sample
sun_samplesheet$diagnosis <- as.factor(sun_samplesheet$diagnosis)

mitogene_dds <- DESeqDataSetFromMatrix(mitogene_counts,
                                       colData = sun_samplesheet,
                                       design = ~ diagnosis)
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
  xlim(-4,4) + ylim(0,22) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_text(nudge_x=.20) +
  scale_colour_identity() +
  geom_vline(xintercept = .6,linetype="dotted") +
  geom_vline(xintercept = -.6,linetype="dotted") +
  labs(title="MitoGene Volcano Plot for IFNa-treated Human T Cells")

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
                                            colData = sun_samplesheet,
                                            design = ~ diagnosis)
mdp_dds_corrected <- DESeq(mdp_dds_corrected)
mdp_res_corrected <-data.frame(results(mdp_dds_corrected))
mdp_res_corrected <- mdp_res_corrected[!is.na(mdp_res_corrected$padj),]
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
  labs(title="MDPSeq Volcano Plot for Leukemia")

mdp_dds_uncorrected <- DESeqDataSetFromMatrix(mdp_counts,
                                              colData = sun_samplesheet,
                                              design = ~ diagnosis)
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
  labs(title="MDPSeq Volcano Plot for Leukemia - Uncorrected")

#let's try this pizza plot thing - partially Brendan's code.
#Pre-processing MDP gtf.
colnames(mdp_gtf) <-
  c("chr",
    "idc",
    "idc2",
    "start",
    "end",
    "idc3",
    "sense",
    "idc4",
    "gene")
mdp_gtf <- separate(mdp_gtf, gene, c("gene", "transcript"), ";")
mdp_gtf <- separate(mdp_gtf, gene, c("trash", "mdp_id"), ">Peptide")
mdp_db <- mdp_gtf[, c("start", "end", "sense", "mdp_id")]
mdp_db <- mdp_db[-c(62, 425), ]

pval <- 0.5
h <- 0.05
borderCol = "gray"

top_results <- mdp_res_corrected[mdp_res_corrected$padj < .05 & abs(mdp_res_corrected$log2FoldChange)>.1,]
top_results$mdp_id <- rownames(top_results) 
list_of_genes <- data.frame(mdp_id = top_results$mdp_id[order(top_results$padj)], 
                          fold = top_results$log2FoldChange[order(top_results$padj)]) %>%
  arrange(fold)
positions <- match(list_of_genes$mdp_id,mdp_db$mdp_id)
positional_table <- merge(drop_na(mdp_db[positions,]),list_of_genes)

positives_fold <- list_of_genes$fold[list_of_genes$fold > 0] %>%
  round(2)
positives_name <- list_of_genes$mdp_id[list_of_genes$fold > 0]

negatives_fold <- list_of_genes$fold[list_of_genes$fold < 0] %>%
  round(2)
negatives_name <- list_of_genes$mdp_id[list_of_genes$fold < 0]

list_of_colors_neg <- colorRampPalette(c("dodgerblue", "yellow"))
list_of_colors_neg <- list_of_colors_neg(length(negatives_name))
list_of_colors_pos <- colorRampPalette(c("yellow", "firebrick"))
list_of_colors_pos <- list_of_colors_pos(length(positives_name))

positives <- data.frame(positives_name,positives_fold,list_of_colors_pos)
negatives <- data.frame(negatives_name,negatives_fold,list_of_colors_neg)

positives_position <- match(positives$positives_name,positional_table$mdp_id)
negatives_position <- match(negatives$negatives_name,positional_table$mdp_id)

positional_table$color <- NULL
positional_table$color[positives_position] <- positives$list_of_colors_pos
positional_table$color[negatives_position] <- negatives$list_of_colors_neg

colors_pos <- data.frame(positives_fold, list_of_colors_pos)
colors_neg <- data.frame(negatives_fold, list_of_colors_neg)
plot(c(-1, 1), c(-1, 1), type = "n", axes = FALSE, ann = FALSE, asp = 1)
draw.sector(0, 360, rou1 = 0.98, rou2 = 0.97, col = "black")

pos <- 0
for (gene in 1:nrow(positional_table)){
  if(gene > 2){
    if(positional_table$start[gene] > positional_table$end[gene-1]){
      pos <- 0
    }  
    if(positional_table$start[gene] < positional_table$end[gene-1]){
      pos <- (pos + 1)
    }
  }
  color <- positional_table$color[gene]
  draw.sector(((positional_table$end[gene]+4300) * 0.021727322), ((positional_table$start[gene]+4300) * 0.021727322),
              col = color, border = borderCol, rou1 = 0.97 - (pos*h), rou2 = (0.97-h) - (pos*h))
}

#These sectors represent the genes
draw.sector(((577 + 4300) * 0.021727322), ((-546 + 4300) * 0.021727322 ), col = "white", border = "black", rou1 = 1.04, rou2 =1)
draw.sector(((647 + 4300) * 0.021727322), ((1601 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((1674 + 4300) * 0.021727322), ((3229 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((3307 + 4300) * 0.021727322), ((4262 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((4470 + 4300) * 0.021727322), ((5511 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((5904 + 4300) * 0.021727322), ((7445 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((7586 + 4300) * 0.021727322), ((8295 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((8364 + 4300) * 0.021727322), ((9207 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((9207 + 4300) * 0.021727322), ((9990 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((10059 + 4300) * 0.021727322), ((10404 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((10470 + 4300) * 0.021727322), ((12138 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((12336 + 4300) * 0.021727322), ((14149 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((14149 + 4300) * 0.021727322), ((14673 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)
draw.sector(((14747 + 4300) * 0.021727322), ((15887 + 4300) * 0.021727322 ), col = "white", clock.wise = FALSE, border = "black", rou1 = 1.04, rou2 = 1)

lgd_ = rep(NA, 11)
lgd_[c(1,6,11)] = c(round(min(negatives_fold),2),0,round(max(positives_fold),2))
legend(x = 1.2, y = 0.5,
       legend = lgd_,
       fill = colorRampPalette(colors = c('dodgerblue','yellow','red'))(11),
       border = NA,
       y.intersp = 0.5,
       cex = 1.5, text.font = 2)

#mess w/ mass spec data
ms_database <- read.csv("~/Desktop/Data/ms_database.csv", comment.char="#")
ms_database$mdp_id <- gsub(">Peptide", "", ms_database$peptide)
joined_database <- merge(ms_database,top_results)
top_hits_sun <- top_results
