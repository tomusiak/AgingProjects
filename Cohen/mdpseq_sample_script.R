#MDPSeq code on constitutively vs. non-constitutively active CARs.
#Based off code provided by Brendan Miller.

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
library(ggrepel)

#Loading in auxiliary code.
setwd("/home/atom/Desktop/AgingProjects/Cohen/")
source("mdpseq_background_correction.R")

setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

#Loading in data.
setwd("/home/atom/Desktop/Data")
list_of_bams <- c("D11OFF.fastq.subread.BAM", "D11ON.fastq.subread.BAM",
                  "D15OFF.fastq.subread.BAM","D15ON.fastq.subread.BAM")
sampleTable<-read.csv("samples.csv")
sampleTable <- sampleTable[3:6,]
bamfiles<-BamFileList(list_of_bams, yieldSize = 2000000)
gtffile <- file.path(getwd(),"mdp_atg_noATC.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
mdp_gtf <-
  read.delim("mdp_atg_noATC.gtf", header = FALSE)
mitogene_gtf <-
  read.delim(
    "mitochondria_db.gtf",
    header = FALSE,
    comment.char = "#"
  )

# Annotating and aligning. Saving & loading optional - summarizeOverlaps() can take time.
#ebg <- exonsBy(txdb, by="gene")
#se <- summarizeOverlaps(features=ebg, reads=list_of_bams,
#                        mode="Union",
#                        singleEnd=TRUE,
#                        ignore.strand=FALSE,
#                        fragments=FALSE,
#                        inter.feature=FALSE) ###this will count each)
#saveRDS(se, "mito_se.rds")
se <- readRDS("mito_se.rds")
mdp_counts <- assays(se)$counts
rownames(mdp_counts) <- sub('.Peptide', '', rownames(mdp_counts))
mdp_counts <- mdp_counts[c(-593,-595),]
keep <- rowSums((mdp_counts)) >= 200 #Removes genes with low counts.
mdp_counts<- mdp_counts[keep,]

#gene_counts <- featureCounts(list_of_bams,annot.ext="~/Desktop/Data/mitochondria_db.gtf",
#                             isGTFAnnotationFile=TRUE,nthreads=5) #Count matrix generation
#raw_count_genes <- gene_counts$counts
#rownames(raw_count_genes) <- c("TF", "RNR1", "TV", "RNR2", "TL1",
#                               "ND1", "TI", "TQ", "TM", "ND2",
#                               "TW", "TA", "TN", "TC","TY",
#                               "CO1", "TS","TD", "CO2","TK",
#                               "ATP8","ATP6","CO3", "TG", "ND3",
#                               "TR","ND4L","ND4","TH","TS2",
#                               "TL2","ND5","ND6","TE","CYB",
#                             "TT","TP")
#write.csv(raw_count_genes,"gene_counts_raw.csv") #Writes count matrix for easier future loading.
raw_count_genes <- read.csv("gene_counts_raw.csv", row.names=1) 
rownames(raw_count_genes) <- sub('MT-', '', rownames(raw_count_genes))
encompass_table <- generateEncompassTable(mdp_gtf,mitogene_gtf)
keep_mitogenes <- c("RNR1","RNR2", "ND1","ND2","CO1","CO2","ATP8",
                    "ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")
encompass_table <- encompass_table[encompass_table$mitogene %in% keep_mitogenes,]
background_table <- determineBackgroundSignal(raw_count_genes)
corrected_counts <- data.matrix(performBackgroundCorrection(background_table,mdp_counts,encompass_table))

####PERFORM DESEQ BETWEEN GROUPS - corrected
col_data <- DataFrame(sampleTable) ###for the diagnosis
col_data$Day <- as.factor(col_data$Day)
col_data$CAR <- as.factor(col_data$CAR)
colnames(corrected_counts) <- col_data$File
colnames(mdp_counts) <- col_data$File
dds_corrected <- DESeqDataSetFromMatrix(corrected_counts,
                              colData = col_data,
                              design = ~ Day + CAR)
dds_corrected <- DESeq(dds_corrected, betaPrior=FALSE)
res_corrected<-results(dds_corrected, name="CAR_ON_vs_OFF")
vsd_corrected <- varianceStabilizingTransformation(dds_corrected) ###mdp-seq script
resOrdered_corrected <- res_corrected[order(abs(res_corrected$pvalue)),]
top_corrected <- head(resOrdered_corrected, 30)
#14A = humanin
#9A = MOTSC
#saveRDS(dds, "shadel_dds.rds")
#dds <- readRDS('shadel_dds.rds')

#write.csv(res, "shadel_results.csv")
#res <- read.csv('shadel_results.csv')

#results <- as.data.frame(res)

#DESeq not corrected
dds <- DESeqDataSetFromMatrix(mdp_counts,
                                        colData = col_data,
                                        design = ~ Day + CAR)
dds <- DESeq(dds)
res <-results(dds, name="CAR_ON_vs_OFF")
vsd <- varianceStabilizingTransformation(dds) ###mdp-seq script
resOrdered <- res[order(abs(res$pvalue)),]
top <- head(resOrdered, 30)
resOrdered["52A",]

#PCa plot not corrected
plotPCA(vsd, "CAR") +
  theme_classic()+
  geom_point(aes(colour = factor(CAR)))+
  xlim(-4,4) +
  ylim(-4,4) +
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

#PCA plot corrected
plotPCA(vsd_corrected, "CAR") +
  
  theme_classic()+
 geom_point(aes(colour = factor(CAR)))+
  xlim(-4,4) +
  ylim(-4,4) +
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

#Count plotting for not-corrected analysis
plotCounts(dds, "32B", intgroup=c("CAR"), returnData=TRUE) %>% 
  ggplot(aes(CAR, count)) + geom_boxplot(aes(color=CAR), lwd = 1.5) + 
  geom_point(aes(color=CAR),size=3, color = "grey") +
  # geom_signif(annotation="pvalue<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_y_log10() + ggtitle("Transcript Differences Between Perma-On CAR and Perma-Control CAR") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Day 11", "Day 15"))+
  theme_minimal()

#Count plotting for corrected analysis
plotCounts(dds_corrected, "32B", intgroup=c("CAR"), returnData=TRUE) %>% 
  ggplot(aes(CAR, count)) + geom_boxplot(aes(color=CAR), lwd = 1.5) + 
  geom_point(aes(color=CAR),size=3, color = "grey") +
  # geom_signif(annotation="pvalue<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_y_log10() + ggtitle("Transcript Differences Between Perma-On CAR and Perma-Control CAR") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Day 11", "Day 15"))+
  theme_minimal()

#Volcano plot for corrected analysis
cols <- densCols(res_corrected$log2FoldChange, -log10(res_corrected$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= res_corrected$log2FoldChange, 
     y = -log10(res_corrected$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=res_corrected$pch, cex=0.4)
gn.selected <- abs(res_corrected$log2FoldChange) >.35 & res_corrected$pvalue < .05
text(res_corrected$log2FoldChange[gn.selected],
     -log10(res_corrected$pvalue)[gn.selected],
     lab=rownames(res_corrected)[gn.selected ], cex=0.6)

#Volcano plot for not corrected analysis.
cols <- densCols(res$log2FoldChange, -log10(res$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))

plot(x= res$log2FoldChange, 
     y = -log10(res$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=res$pch, cex=0.4)
gn.selected <- abs(res$log2FoldChange) >.35 & res$pvalue < .05
text(res$log2FoldChange[gn.selected],
     -log10(res$pvalue)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.6)

raw_count_genes <-raw_count_genes[rownames(raw_count_genes) %in% keep_mitogenes,]
mito_dds <- DESeqDataSetFromMatrix(raw_count_genes,
                                        colData = col_data,
                                        design = ~ Day + CAR)
mito_dds <- DESeq(mito_dds)
mito_res <-results(mito_dds, name="CAR_ON_vs_OFF")
mito_resLFC <- lfcShrink(mito_dds, "CAR_ON_vs_OFF")
mito_vsd <- varianceStabilizingTransformation(mito_dds) ###mdp-seq script
mito_resOrdered <- mito_res[order(abs(mito_res$pvalue)),]
mito_top <- head(mito_resOrdered, 30)

coloring <- encompass_table %>% group_by (mitogene)

#Volcano plot for not corrected analysis.
cols <- densCols(mito_res$log2FoldChange, -log10(mito_res$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= mito_res$log2FoldChange, 
     y = -log10(mito_res$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-2,2),
     ylim=c(0,100),
     pch=mito_resLFC$pch, cex=0.4)
gn.selected <- abs(mito_resLFC$log2FoldChange) >.2 & mito_res$pvalue < .05
text(mito_res$log2FoldChange[gn.selected],
     -log10(mito_res$pvalue)[gn.selected],
     lab=rownames(mito_res)[gn.selected ], cex=0.6)

corrected_counts["241C",]
mdp_counts["241C",]
raw_count_genes["CO1",]
res_corrected["241C",]
res["241C",]
mito_resOrdered["CO1",]

corrected_counts["32B",]
mdp_counts["32B",]
raw_count_genes["CO2",]
res_corrected["32B",]
res["32B",]
mito_res["CO2",]

corrected_counts["106D",]
mdp_counts["106D",]
raw_count_genes["CO3",]
res_corrected["106D",]
res["106D",]
mito_res["CO3",]

#Time to make some plots for a Cohen lab presentation.
cols <- densCols(res$log2FoldChange, -log10(res$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= res$log2FoldChange, 
     y = -log10(res$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot of Exhausted CAR-Ts vs. Control CAR-Ts", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=res$pch, cex=1)
gn.selected <- rownames(res)=="241C" | rownames(res)=="32B"
text(res$log2FoldChange[gn.selected],
     -log10(res$pvalue)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=.8)

co2_counts <- data.frame(raw_count_genes["CO2",] %>% t())
co2_counts$status <- c("Control","Exhausted","Control","Exhausted")
peptide_counts <- data.frame(mdp_counts["32B",])
peptide_counts$status <- c("Control","Exhausted","Control","Exhausted")
peptide_counts_corrected <- data.frame(corrected_counts["32B",])
peptide_counts_corrected$status <-c("Control","Exhausted","Control","Exhausted")
co2_counts$CO2 <- co2_counts$CO2/mean(co2_counts$CO2)
colnames(peptide_counts) <- c("counts","status")
colnames(peptide_counts_corrected) <- c("counts","status")
peptide_counts$counts <- peptide_counts$counts/mean(peptide_counts$counts)
peptide_counts_corrected$counts <- peptide_counts_corrected$counts/mean(peptide_counts_corrected$counts)
co2_summary <- summarySE(co2_counts,"CO2", "status")
peptide_counts_summary <- summarySE(peptide_counts,"counts", "status")
peptide_counts_corrected_summary<- summarySE(peptide_counts_corrected,"counts", "status")

ggplot(co2_summary, aes(x=status, y=CO2)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill="pink") +
  ylim(0,2) +
  geom_errorbar(aes(ymin=CO2-se,ymax=CO2+se),position = position_dodge(width = .9),width=.2) + 
  theme_classic() + 
  labs(title="CO2 Expression on Control vs. Exhausted T Cells",x = "Control or Exhausted",
       y="Normalized Counts") +
  scale_fill_manual(values=c("pink"))

ggplot(peptide_counts_corrected_summary, aes(x=status, y=counts)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill="magenta") +
  ylim(0,2) +
  geom_errorbar(aes(ymin=counts-se,ymax=counts+se),position = position_dodge(width = .9),width=.2) + 
  theme_classic() + 
  labs(title="32B Expression on Control vs. Exhausted T Cells",x = "Control or Exhausted",
       y="Normalized Counts") +
  scale_fill_manual(values=c("pink"))

co1_counts <- data.frame(raw_count_genes["CO1",] %>% t())
co1_counts$status <- c("Control","Exhausted","Control","Exhausted")
peptide_counts <- data.frame(mdp_counts["241C",])
peptide_counts$status <- c("Control","Exhausted","Control","Exhausted")
peptide_counts_corrected <- data.frame(corrected_counts["241C",])
peptide_counts_corrected$status <-c("Control","Exhausted","Control","Exhausted")
co1_counts$CO1 <- co1_counts$CO1/mean(co1_counts$CO1)
colnames(peptide_counts) <- c("counts","status")
colnames(peptide_counts_corrected) <- c("counts","status")
peptide_counts$counts <- peptide_counts$counts/mean(peptide_counts$counts)
peptide_counts_corrected$counts <- peptide_counts_corrected$counts/mean(peptide_counts_corrected$counts)
co1_summary <- summarySE(co1_counts,"CO1", "status")
peptide_counts_summary <- summarySE(peptide_counts,"counts", "status")
peptide_counts_corrected_summary<- summarySE(peptide_counts_corrected,"counts", "status")

ggplot(co1_summary, aes(x=status, y=CO1)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill="pink") +
  ylim(0,2) +
  geom_errorbar(aes(ymin=CO1-se,ymax=CO1+se),position = position_dodge(width = .9),width=.2) + 
  theme_classic() + 
  labs(title="CO1 Expression on Control vs. Exhausted T Cells",x = "Control or Exhausted",
       y="Normalized Counts") +
  scale_fill_manual(values=c("red"))

ggplot(peptide_counts_corrected_summary, aes(x=status, y=counts)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), fill="magenta") +
  ylim(0,2) +
  geom_errorbar(aes(ymin=counts-se,ymax=counts+se),position = position_dodge(width = .9),width=.2) + 
  theme_classic() + 
  labs(title="241C Expression on Control vs. Exhausted T Cells",x = "Control or Exhausted",
       y="Normalized Counts") +
  scale_fill_manual(values=c("red"))

cols <- densCols(res_corrected$log2FoldChange, -log10(res_corrected$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= res_corrected$log2FoldChange, 
     y = -log10(res_corrected$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot of Exhausted CAR-Ts vs. Control CAR-Ts", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=res_corrected$pch, cex=1)
gn.selected <- rownames(res_corrected)=="241C" | rownames(res_corrected)=="32B"
text(res_corrected$log2FoldChange[gn.selected],
     -log10(res_corrected$pvalue)[gn.selected],
     lab=rownames(res_corrected)[gn.selected ], cex=.8)
