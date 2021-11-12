##########====HEPATIC SHY5Y MDP-SEQ SCRIPT=====#########
setwd("/home/atom/Desktop/Data/")
list_of_bams <- c("D7OFF.fastq.subread.BAM","D7ON.fastq.subread.BAM","D11OFF.fastq.subread.BAM",
                  "D11ON.fastq.subread.BAM","D15OFF.fastq.subread.BAM","D15ON.fastq.subread.BAM")
sampleTable<-read.csv("samples.csv")

slist.files(getwd())
filenames <- file.path(getwd(), paste0(sampleTable$id, "Aligned.out_sorted_chrM.bam"))
file.exists(filenames)
#####LOAD BAM FILES----
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
bamfiles<-BamFileList(list_of_bams, yieldSize = 2000000)
seqinfo(bamfiles[1])

###LOAD
library("Biobase")
library("AnnotationDbi")
library("GenomicFeatures")
gtffile <- file.path(getwd(),"mdp_atg_noATC.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

######annotate exons (or MDPs)
ebg <- exonsBy(txdb, by="gene")
ebg

#####ALIGN TO MDP CALLS
library("matrixStats")
library("DelayedArray")
library("SummarizedExperiment")
library("GenomicAlignments")
library("BiocParallel")

se <- summarizeOverlaps(features=ebg, reads=list_of_bams,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=FALSE,
                        fragments=FALSE,
                        inter.feature=FALSE) ###this will count each)

#saveRDS(se, "shadel_se_all.rds")
se <- readRDS('shadel_se_all.rds')

####PERFORM DESEQ BETWEEN GROUPS
library("DESeq2")
library(dplyr)
colData(se) <- DataFrame(sampleTable) ###for the diagnosis
colData(se)$Day <- as.factor(colData(se)$Day)
dds <- DESeqDataSet(se, design = ~  CAR + Day)
dds <- DESeq(dds)
#saveRDS(dds, "shadel_dds.rds")
dds <- readRDS('shadel_dds.rds')

res <- results(dds)
#write.csv(res, "shadel_results.csv")
#res <- read.csv('shadel_results.csv')

#results <- as.data.frame(res)

# #######GENERATEA FIGURES
library("pheatmap")
library(ggplot2)
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
plotPCA(vsd, "Day") +
  theme_classic()+
 geom_point(aes(colour = factor(Day)))+
  geom_label_repel(aes(label = Day),
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
resOrdered <- res[order(res$pvalue),]
top <- head(resOrdered, 30)
top[">Peptide3A",]
topGene <- ">Peptide98D"
#14A = humanin
#9A = MOTSC

plotCounts(dds, topGene, intgroup=c("Day","CAR"), returnData=TRUE) %>% 
  ggplot(aes(Day, count)) + geom_boxplot(aes(color=CAR), lwd = 1.5) + 
  geom_point(aes(color=CAR),size=3, color = "grey") +
  # geom_signif(annotation="pAdj<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_fill_manual(values=colours, name="Group")+
  scale_color_manual(values=colours)+
  scale_y_log10() + ggtitle("Transcript Differences Between Perma-On CAR and Perma-Off CAR") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Day 7","Day 11", "Day 15"))+
  theme_minimal()


