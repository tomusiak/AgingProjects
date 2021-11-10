##########====HEPATIC SHY5Y MDP-SEQ SCRIPT=====#########
setwd("/home/atom/Desktop/Data/")
list_of_bams <- c("SRR13961047.1.fastq.subread.BAM","SRR13961048.1.fastq.subread.BAM","SRR13961049.1.fastq.subread.BAM","SRR13961050.1.fastq.subread.BAM")
sampleTable<-read.csv("shadel_sample_table.csv")

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
colData(se)
se$disease <- relevel(se$passage, "Passage5")
dds <- DESeqDataSet(se, design = ~  passage)
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
mat <- assay(vsd)[ head(order(res$padj),100), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[,c("id","passage")])
pheatmap(mat, annotation_col=df)

library(ggrepel)
plotPCA(vsd, "id") +
  theme_classic()+
 geom_point(aes(colour = factor(id)))+
  geom_label_repel(aes(label = id),
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




library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                
                xlim=c(-6,6),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Healthy vs. NASH',
                titleLabSize = 50,
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.1,
                FCcutoff = 0.5,
                transcriptPointSize = 5.5,
                transcriptLabSize = 6.0,
                col=c('black', 'red', 'orange', 'blue'),
                legendPosition = 'right',
                legend=c('Not Significant','Not Significant, But Big Fold Change',
                         'Signficant, Minimal Fold Change','Signficant, Big Fold Change'),
                legendLabSize = 16,
                legendIconSize = 10.0,    transcriptLabCol = 'black',
                transcriptLabFace = 'bold') +
  labs(title = "Senescence MDP-Seq Fold Changes",
       subtitle = "Each dot represents one smORF",
       caption = "", 
       x = "Log2 Fold Change", y = "-Log10 adjusted P")

topGene <- rownames(res)[148] #humanin
topGene <- rownames(res)[644] #shlp2
topGene <- rownames(res)[525] #shmoose
plotCounts(dds, topGene, "passage", xlab="Group", ylabel="Normalized Counts") 

colours <- c("blue", "red")
plotCounts(dds, topGene, intgroup="passage", returnData=TRUE) %>% 
  ggplot(aes(passage, count)) + geom_boxplot(aes(fill=passage), lwd = 1.5) + 
  geom_point(aes(fill=passage),size=3, color = "grey") +
  # geom_signif(annotation="pAdj<0.15", textsize = 5.5,
  #             y_position=4.76, xmin=1, xmax=4, tip_length = 0.1, size = 1, fontface="bold", color="black") +
  scale_fill_manual(values=colours, name="Group")+
  scale_color_manual(values=colours)+
  scale_y_log10() + ggtitle("SHMOOSE Transcript Differences Between Senescence") + 
  ylab("Normalized Count")+
  scale_x_discrete(labels= c("Control-SHMOOSE WT (n=59)","AD-SHMOOSE WT (n=59)", "Control-SHMOOSE WT (n=59)"))+
  theme_minimal()+
  theme(axis.text=element_text(size=30),
        plot.title=element_text(size=35, face="bold"),
        axis.title=element_text(size=20,face="bold"), 
        legend.text=element_text(size=18),
        legend.title=element_text(size=20, face="bold"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=20, face="bold", color="black"))

