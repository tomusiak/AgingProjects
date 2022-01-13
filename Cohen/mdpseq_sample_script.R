#MDPSeq code on constitutively vs. non-constitutively active CARs.
#Based off code provided by Brendan Miller.

#Loading in background correction code.
source("mdpseq_background_correction.R")

#Loading in libraries.
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

#Loading in data.
setwd("/home/atom/Desktop/Data/")
list_of_bams <- c("D11OFF.fastq.subread.BAM", "D11ON.fastq.subread.BAM",
                  "D15OFF.fastq.subread.BAM","D15ON.fastq.subread.BAM")
sampleTable<-read.csv("samples.csv")
sampleTable <- sampleTable[3:6,]
bamfiles<-BamFileList(list_of_bams, yieldSize = 2000000)
gtffile <- file.path(getwd(),"mdp_atg_noATC.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())

# Annotating and aligning. Saving & loading optional - summarizeOverlaps() can take time.
ebg <- exonsBy(txdb, by="gene")
#se <- summarizeOverlaps(features=ebg, reads=list_of_bams,
#                        mode="Union",
#                        singleEnd=TRUE,
#                        ignore.strand=FALSE,
#                        fragments=FALSE,
#                        inter.feature=FALSE) ###this will count each)
#saveRDS(se, "mito_se.rds")
se <- readRDS('mito_se.rds')
mdp_counts <- assays(se)$counts

#gene_counts <- featureCounts(list_of_bams,annot.ext="~/Desktop/Data/mitochondria_db.gtf",
#                             isGTFAnnotationFile=TRUE,nthreads=5) #Count matrix generation
#raw_count_genes <- gene_counts$counts
#rownames(raw_count_genes) <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1",
#                               "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2",
#                               "MT-TW", "MT-TA", "MT-TN", "MT-TC","MT-TY",
#                               "MT-CO1", "MT-TS","MT-TD", "MT-CO2","MT-TK",
#                               "MT-ATP8","MT-ATP6","MT-CO3", "MT-TG", "MT-ND3",
#                               "MT-TR","MT-ND4L","MT-ND4","MT-TH","MT-TS2",
#                               "MT-TL2","MT-ND5","MT-ND6","MT-TE","MT-CYB",
#                             "MT-TT","MT-TP")
#write.csv(raw_count_genes,"gene_counts_raw.csv") #Writes count matrix for easier future loading.
raw_count_genes <- read.csv("gene_counts_raw.csv", row.names=1) 


####PERFORM DESEQ BETWEEN GROUPS
colData(se) <- DataFrame(sampleTable) ###for the diagnosis
colData(se)$Day <- as.factor(colData(se)$Day)
dds <- DESeqDataSet(se, design = ~  CAR + Day)
dds <- DESeq(dds)
#saveRDS(dds, "shadel_dds.rds")
#dds <- readRDS('shadel_dds.rds')

res <- results(dds)
#write.csv(res, "shadel_results.csv")
#res <- read.csv('shadel_results.csv')

#results <- as.data.frame(res)

# #######GENERATEA FIGURES
vsd <- varianceStabilizingTransformation(dds) ###mdp-seq script

mat <- assay(vsd)
on_minus_off_day7 <- mat[,2] - mat[,1]
on_minus_off_day11 <- mat[,4] - mat[,3]
on_minus_off_day15 <- mat[,6] - mat[,5]
data <- data.frame(on_minus_off_day7)
data$day11 <- on_minus_off_day11
data$day15 <- on_minus_off_day15
df <- as.data.frame(colData(vsd)[,c("Day","CAR")])
rownames(new_annot) <- as.data.frame(df[c(1,3,5),1])
rownames(new_annot) <- c(7,11,15)
colnames(data) <- c(7,11,15)
pheatmap(data, annotation_col=new_annot,axis=1)

data
processing <- as.data.frame((data$"15") + (data$"11") + (data$"7")) 
rownames(processing) <- rownames(data)
colnames(processing) <- "values"
head(processing)
sorted <- processing[order(processing$"values"),,drop=FALSE]
head(sorted,10)
tail(sorted,20)

library(ggrepel)
plotPCA(vsd, "CAR") +
  theme_classic()+
 geom_point(aes(colour = factor(CAR)))+
  geom_label_repel(aes(label = CAR),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  labs(title = "Senescence MDP-Seq ",
       subtitle = "Each dot represents one donor",
       #caption = "Date: 03-21-2019; P17", 
       x = "PC1 (50%)", y = "PC2(21%)") +
  theme(axis.text=element_text(size=30),
        plot.title=element_text(size=35, face="bold"),
        plot.subtitle=element_text(size=20, hjust=0.0, face="italic", color="black"),
        plot.caption=element_text(size=14, face="italic", color="black"),
        axis.title=element_text(size=25,face="bold"), 
        legend.text=element_text(size=25),
        legend.title=element_text(size=25, face="bold"),
        #legend.key.size = unit(2.5, 'lines'),
        axis.title.x=element_text(size=25, face="bold", color="black"),
        axis.text.x=element_text(size=15, face="bold", color="black"),
        axis.text.y=element_text(size=15, face="bold", colour = "black"))

res<-results(dds, name="CAR_ON_vs_OFF")
res<-results(dds,contrast=c("CAR","ON","OFF"))
resLFC <- lfcShrink(dds, coef="CAR_ON_vs_OFF", type="apeglm")
resOrdered <- res[order(res$log2FoldChange),]
top <- head(resOrdered, 30)
gene <- resOrdered[">Peptide9A",]
#14A = humanin
#9A = MOTSC

plotCounts(dds, ">Peptide9A", intgroup=c("CAR"), returnData=TRUE) %>% 
  ggplot(aes(CAR, count)) + geom_boxplot(aes(color=CAR), lwd = 1.5) + 
  geom_point(aes(color=CAR),size=3, color = "grey") +
  # geom_signif(annotation="pAdj<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_y_log10() + ggtitle("Transcript Differences Between Perma-On CAR and Perma-Off CAR") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Day 11", "Day 15"))+
  theme_minimal()

cols <- densCols(resLFC$log2FoldChange, -log10(resLFC$padj),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= resLFC$log2FoldChange, 
     y = -log10(resLFC$padj), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-1,1),
     ylim=c(0,4),
     pch=resLFC$pch, cex=0.4)
gn.selected <- abs(resLFC$log2FoldChange) >.35 & resLFC$padj < .05
text(resLFC$log2FoldChange[gn.selected],
     -log10(resLFC$padj)[gn.selected],
     lab=rownames(resLFC)[gn.selected ], cex=0.6)
