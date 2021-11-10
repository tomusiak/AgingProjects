library(Rsamtools)
library(Rsubread)
library(DESeq2)
library(apeglm)
library(magrittr)
setwd("/home/atom/Desktop/Data/")
#buildindex(basename="reference_index_mouse",reference="GRCm38.primary_assembly.genome.fa",memory=6000)
#align("reference_index_mouse", "SRR9732704.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732705.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732706.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732707.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732708.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732709.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732710.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732711.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732712.1.fastq", nthreads=5)
#align("reference_index_mouse", "SRR9732713.1.fastq", nthreads=5)
list_of_sams <- c("SRR9732704.1.fastq.subread.BAM","SRR9732705.1.fastq.subread.BAM","SRR9732706.1.fastq.subread.BAM","SRR9732707.1.fastq.subread.BAM",
                "SRR9732708.1.fastq.subread.BAM", "SRR9732709.1.fastq.subread.BAM","SRR9732710.1.fastq.subread.BAM","SRR9732711.1.fastq.subread.BAM",
                "SRR9732712.1.fastq.subread.BAM","SRR9732713.1.fastq.subread.BAM")
counts <- featureCounts(list_of_sams,annot.inbuilt="mm10",useMetaFeatures=TRUE)
orig_matrix <- counts$counts
colnames(orig_matrix) <- c("CTRL_1","CTRL_2","CTRL_3","CTRL_4","CTRL_5",
                           "MOTSC_1","MOTSC_2","MOTSC_3","MOTSC_4","MOTSC_5")
orig_data_metadata <- read.csv("metadata - Sheet1.csv")
rownames(orig_data_metadata) <- orig_data_metadata$X
head(orig_matrix)
de_matrix <- DESeqDataSetFromMatrix(countData = orig_matrix,
                              colData = orig_data_metadata,
                              design = ~ condition)
de_matrix <- DESeq(de_matrix)
res<-results(de_matrix, name="condition_MOTSC_vs_CTRL")
res<-results(de_matrix,contrast=c("condition","MOTSC","CTRL"))
resLFC <- lfcShrink(de_matrix, coef="condition_MOTSC_vs_CTRL", type="apeglm")
resOrdered <- res[order(res$pvalue),]
top <- head(resOrdered, 80)
top["16476",]
plotMA(resLFC, ylim=c(-2,2))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
plotCounts(de_matrix, gene="12675", intgroup="condition")
d <- plotCounts(de_matrix, gene="12675", intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + theme_classic() + 
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd, width=.2)) +
  scale_y_continuous(name="Normalized Counts", limits = c(100,200))
d$mean <- c(121,121,121,121,121,157,157,157,157,157)
d$sd <- c(6,6,6,6,6,14,14,14,14,14)

#buildindex(basename="reference_index",reference="GRCh38.primary_assembly.genome.fa",memory=6000,gappedIndex=TRUE)
#align("reference_index", "SRR13961047.1.fastq", nthreads=5)
#align("reference_index", "SRR13961048.1.fastq", nthreads=5)
#align("reference_index", "SRR13961049.1.fastq", nthreads=5)
#align("reference_index", "SRR13961050.1.fastq", nthreads=5)
list_of_sams_mdpseq <- c("SRR13961047.1.fastq.subread.BAM","SRR13961048.1.fastq.subread.BAM",
                         "SRR13961049.1.fastq.subread.BAM","SRR13961050.1.fastq.subread.BAM")
mdpseq_counts <- featureCounts(list_of_sams_mdpseq, annot.ext="mdp_atg_noATC.gtf",isGTFAnnotationFile=TRUE)
