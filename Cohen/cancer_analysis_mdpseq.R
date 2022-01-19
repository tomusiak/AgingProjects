keep_mitogenes <- c("MT-RNR1","MT-RNR2", "MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8",
                    "MT-ATP6","MT-CO3","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MY-CYB")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE, removeOutliers=F, logt=F) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         
                         if (removeOutliers) {
                           print('pre')
                           print(xx[[col]])
                           #xx[[col]] <- boxB(xx[[col]])
                           print('post')
                           #print(exp(rm.outlier(log(xx[[col]]), fill = FALSE, median = FALSE, opposite = FALSE)))
                           print(scores(log(xx[[col]]), type="t", prob=0.90) )
                           #print(xx[[col]][boxB(xx[[col]], method="resistant", logt=logt, k=2.5)[["outliers"]]])
                         }
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

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

setwd("/home/atom/Desktop/Data/cancer_data")
#align(index="../reference_index",readfile1="SRR4427174.fastq",output_file="SRR4427174.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427175.fastq",output_file="SRR4427175.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427176.fastq",output_file="SRR4427176.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427177.fastq",output_file="SRR4427177.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427178.fastq",output_file="SRR4427178.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427179.fastq",output_file="SRR4427179.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427180.fastq",output_file="SRR4427180.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427181.fastq",output_file="SRR4427181.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427182.fastq",output_file="SRR4427182.bam",nthreads=4)
#align(index="../reference_index",readfile1="SRR4427183.fastq",output_file="SRR4427183.bam",nthreads=4)
sample_table<-read.csv("patient_data.csv")[1:8,]
bams <- list.files(pattern = "*.bam$")

setwd("..")
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

setwd("../Data/cancer_data")
ebg <- exonsBy(txdb, by="gene")
#se_cancer <- summarizeOverlaps(features=ebg, reads=bams,
#                        mode="Union",
#                        singleEnd=TRUE,
#                        ignore.strand=FALSE,
#                        fragments=FALSE,
#                        inter.feature=FALSE) ###this will count each)
#saveRDS(se_cancer, "mdp_se_cancer.rds")
se_cancer <- readRDS("mdp_se_cancer.rds")
mdp_counts <- assays(se_cancer)$counts
rownames(mdp_counts) <- sub('.Peptide', '', rownames(mdp_counts))
mdp_counts <- mdp_counts[c(-593,-595),]
keep <- rowSums((mdp_counts)) >= 200 #Removes genes with low counts.
mdp_counts<- mdp_counts[keep,1:8]

#mitogene_counts <- featureCounts(bams,annot.ext="~/Desktop/Data/mitochondria_db.gtf",
#                             isGTFAnnotationFile=TRUE,nthreads=5) #Count matrix generation
#mitogene_count_matrix <- mitogene_counts$counts
#rownames(mitogene_count_matrix) <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1",
#                               "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2",
#                               "MT-TW", "MT-TA", "MT-TN", "MT-TC","MT-TY",
#                               "MT-CO1", "MT-TS","MT-TD", "MT-CO2","MT-TK",
#                               "MT-ATP8","MT-ATP6","MT-CO3", "MT-TG", "MT-ND3",
#                               "MT-TR","MT-ND4L","MT-ND4","MT-TH","MT-TS2",
#                               "MT-TL2","MT-ND5","MT-ND6","MT-TE","MT-CYB",
#                             "MT-TT","MT-TP")
#write.csv(mitogene_count_matrix,"mitogene_count_matrix.csv") #Writes count matrix for easier future loading.
mitogene_count_matrix <- read.csv("mitogene_count_matrix.csv", row.names=1)[,1:8]

setwd("/home/atom/Desktop/AgingProjects/Cohen/")
source("mdpseq_background_correction.R")

mitogene_count_matrix <-mitogene_count_matrix[rownames(mitogene_count_matrix) %in% keep_mitogenes,]
mitogene_dds <- DESeqDataSetFromMatrix(mitogene_count_matrix,
                                   colData = sample_table,
                                   design = ~ donor + status)
mitogene_dds <- DESeq(mitogene_dds)
mitogene_res <-results(mitogene_dds, name="status_tumor_vs_normal")
mitogene_vsd <- varianceStabilizingTransformation(mitogene_dds) ###mdp-seq script
mitogene_resOrdered <- mitogene_res[order(abs(mitogene_res$pvalue)),]
mitogene_top <- head(mitogene_resOrdered, 30)

cols <- densCols(mitogene_res$log2FoldChange, -log10(mitogene_res$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= mitogene_res$log2FoldChange, 
     y = -log10(mitogene_res$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=mitogene_res$pch, cex=0.4)
gn.selected <- abs(mitogene_res$log2FoldChange) >.2 & mitogene_res$pvalue < .05
text(mitogene_res$log2FoldChange[gn.selected],
     -log10(mitogene_res$pvalue)[gn.selected],
     lab=rownames(mitogene_res)[gn.selected ], cex=0.6)

encompass_table <- generateEncompassTable(mdp_gtf,mitogene_gtf)
encompass_table <- encompass_table[encompass_table$mitogene %in% keep_mitogenes,]
background_table <- determineBackgroundSignal(mitogene_count_matrix)
corrected_counts <- data.matrix(performBackgroundCorrection(background_table,mdp_counts,encompass_table))

#not corrected
mdp_dds <- DESeqDataSetFromMatrix(mdp_counts,
                              colData = sample_table,
                              design = ~ donor + status)
mdp_dds <- DESeq(mdp_dds)
mdp_res <-results(mdp_dds, name="status_tumor_vs_normal")
mdp_vsd <- varianceStabilizingTransformation(mdp_dds) ###mdp-seq script
mdp_resOrdered <- mdp_res[order(abs(mdp_res$pvalue)),]
mdp_top <- head(mdp_resOrdered, 30)
mdp_resOrdered["173C",]

#corrected
mdp_dds_corrected <- DESeqDataSetFromMatrix(corrected_counts,
                                  colData = sample_table,
                                  design = ~ donor + status)
mdp_dds_corrected <- DESeq(mdp_dds_corrected)
mdp_res_corrected <-results(mdp_dds_corrected, name="status_tumor_vs_normal")
mdp_vsd_corrected <- varianceStabilizingTransformation(mdp_dds_corrected) ###mdp-seq script
mdp_resOrdered_corrected <- mdp_res_corrected[order(abs(mdp_res_corrected$pvalue)),]
mdp_top_corrected <- head(mdp_resOrdered_corrected, 30)
mdp_resOrdered_corrected["192C",]

plotPCA(mdp_vsd, "status") +
  theme_classic()+
  geom_point(aes(colour = factor(status)))+
  xlim(-15,15) +
  ylim(-15,15) +
  geom_label_repel(aes(label = status),
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

cols <- densCols(mdp_res_corrected$log2FoldChange, -log10(mdp_res_corrected$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= mdp_res_corrected$log2FoldChange, 
     y = -log10(mdp_res_corrected$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=mdp_res_corrected$pch, cex=0.4)
gn.selected <- abs(mdp_res_corrected$log2FoldChange) >.5 & mdp_res_corrected$pvalue < .05
text(mdp_res_corrected$log2FoldChange[gn.selected],
     -log10(mdp_res_corrected$pvalue)[gn.selected],
     lab=rownames(mdp_res_corrected)[gn.selected ], cex=0.6)

cols <- densCols(mdp_res$log2FoldChange, -log10(mdp_res$pvalue),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= mdp_res$log2FoldChange, 
     y = -log10(mdp_res$pvalue), 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     xlim=c(-2,2),
     ylim=c(0,10),
     pch=mdp_res$pch, cex=0.4)
gn.selected <- abs(mdp_res$log2FoldChange) >.5 & mdp_res$pvalue < .05
text(mdp_res$log2FoldChange[gn.selected],
     -log10(mdp_res$pvalue)[gn.selected],
     lab=rownames(mdp_res)[gn.selected ], cex=0.6)

#playtime
head(mdp_resOrdered_corrected,15)
corrected_counts["53D",]
mdp_counts["53D",]
mdp_res_corrected["53D",]
mdp_res["53D",]
mitogene_resOrdered["MT-RNR2",]
rnr2_counts <- data.frame(mitogene_count_matrix["MT-RNR2",] %>% t())
rnr2_counts$status <- c("normal","tumor","normal","tumor","normal","tumor","normal","tumor")
peptide_counts <- data.frame(mdp_counts["53D",])
peptide_counts$status <- c("normal","tumor","normal","tumor","normal","tumor","normal","tumor")
peptide_counts_corrected <- data.frame(corrected_counts["53D",])
peptide_counts_corrected$status <- c("normal","tumor","normal","tumor","normal","tumor","normal","tumor")
rnr2_counts$MT.RNR2 <- rnr2_counts$MT.RNR2/mean(rnr2_counts$MT.RNR2)
colnames(peptide_counts) <- c("counts","status")
colnames(peptide_counts_corrected) <- c("counts","status")
peptide_counts$counts <- peptide_counts$counts/mean(peptide_counts$counts)
peptide_counts_corrected$counts <- peptide_counts_corrected$counts/mean(peptide_counts_corrected$counts)
rnr2_summary <- summarySE(rnr2_counts,"MT.RNR2", "status")
peptide_counts_summary <- summarySE(peptide_counts,"counts", "status")
peptide_counts_corrected_summary<- summarySE(peptide_counts_corrected,"counts", "status")

ggplot(rnr2_summary, aes(x=status, y=MT.RNR2)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  ylim(0,1.5) +
  geom_errorbar(aes(ymin=MT.RNR2-se,ymax=MT.RNR2+se),position = position_dodge(width = .9),width=.2) + 
  theme_classic() + 
  labs(title="KLRG1 Expression On T Cell Subsets in SA-bGal High Or Low CD8+ CM Cells",x = "Younger or Older",
       y="Normalized Counts",
       fill="SAbGal High or Low") +
  scale_fill_manual(values=c("pink","firebrick4"))

ggplot(peptide_counts_summary, aes(x=status, y=counts)) +
  ylim(0,1.5) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=counts-se,ymax=counts+se),position = position_dodge(width = .9),width=.2) + 
  theme_classic() + 
  labs(title="KLRG1 Expression On T Cell Subsets in SA-bGal High Or Low CD8+ CM Cells",x = "Younger or Older",
       y="Normalized Counts",
       fill="SAbGal High or Low") +
  scale_fill_manual(values=c("pink","firebrick4"))

ggplot(peptide_counts_corrected_summary, aes(x=status, y=counts)) +
  ylim(0,1.5) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=counts-se,ymax=counts+se),position = position_dodge(width = .9),width=.2) + 
  theme_classic() + 
  labs(title="KLRG1 Expression On T Cell Subsets in SA-bGal High Or Low CD8+ CM Cells",x = "Younger or Older",
       y="Normalized Counts",
       fill="SAbGal High or Low") +
  scale_fill_manual(values=c("pink","firebrick4"))

plotCounts(mitogene_dds, "MT-RNR2", intgroup=c("status"), returnData=TRUE) %>% 
  ggplot(aes(status, count)) + geom_boxplot(aes(color=status), lwd = 1.5) + 
  geom_point(aes(color=status),size=3, color = "grey") +
  # geom_signif(annotation="pvalue<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_y_log10() + ggtitle("Transcript Differences Between Perma-On CAR and Perma-Off CAR") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Day 11", "Day 15"))+
  theme_minimal()

plotCounts(mdp_dds, "53D", intgroup=c("status"), returnData=TRUE) %>% 
  ggplot(aes(status, count)) + geom_boxplot(aes(color=status), lwd = 1.5) + 
  geom_point(aes(color=status),size=3, color = "grey") +
  # geom_signif(annotation="pvalue<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_y_log10() + ggtitle("Transcript Differences Between Perma-On CAR and Perma-Off CAR") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Day 11", "Day 15"))+
  theme_minimal()

plotCounts(mdp_dds_corrected, "53D", intgroup=c("status"), returnData=TRUE) %>% 
  ggplot(aes(status, count)) + geom_boxplot(aes(color=status), lwd = 1.5) + 
  geom_point(aes(color=status),size=3, color = "grey") +
  # geom_signif(annotation="pvalue<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_y_log10() + ggtitle("Transcript Differences Between Perma-On CAR and Perma-Off CAR") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Day 11", "Day 15"))+
  theme_minimal()

corrected_counts["176C",]
mdp_counts["176C",]
mdp_res_corrected["176C",]
mdp_res["176C",]
mitogene_resOrdered["MT-RNR2",]

corrected_counts["51D",]
mdp_counts["51D",]
mitogene_count_matrix["MT-RNR2",]
mdp_res_corrected["51D",]
mdp_res["51D",]
mitogene_resOrdered["MT-RNR2",]

corrected_counts["36A",]
mdp_counts["36A",]
mitogene_count_matrix["MT-CO1",]
mdp_res_corrected["36A",]
mdp_res["36A",]
mitogene_resOrdered["MT-CO1",]

corrected_counts["53D",]
mdp_counts["53D",]
mitogene_count_matrix["MT-RNR2",]
mdp_res_corrected["53D",]
mdp_res["53D",]
mitogene_resOrdered["MT-RNR2",]
