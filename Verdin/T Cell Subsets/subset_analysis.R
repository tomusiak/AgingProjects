#Grabs some useful scripts.
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

#Sets location of data 
setwd("/home/atom/Desktop/Data/subset_data") #Sets directory.

#Libraries to import.
library(WGCNA)
library(missMethyl)
library(limma)
library(dplyr)
require("biomaRt")
library("RColorBrewer")
library(DMRcate)
library(readxl)
library(ggplot2)
library(GOfuncR)
library(readr)
library(plotrix)
library(sesame)
library(umap)
library(splines)
library(gam)
require(clusterExperiment)
library(stringr)
library(tidyr)


clock_data <- read.csv("~/Desktop/Data/subset_data/clock_data.csv")
subset_metadata <- read.csv("~/Desktop/Data/subset_data/subset_metadata.csv")

#D3 was obviously mishandled, removing from dataset...
clock_data <- clock_data[clock_data$SampleID!="D3"]
subset_metadata <- subset_metadata[subset_metadata$SampleID != "D3",]

all_data <- merge(clock_data,subset_metadata)
all_data$type <- as.character(all_data$type)
all_data$type <- factor(all_data$type, levels=c("naive", "central_memory", "effector_memory","temra"))
subset_data <- all_data[all_data$sabgal_sample==FALSE,]
subset_data$diff <- subset_data$DNAmAge-subset_data$age

summary <- getSummary(subset_data,"diff", "type")
summary$type <- c("Naive","Central Memory","Effector Memory","TEMRA")
summary$type <- factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=all_data, aes(x=as.factor(final_dna), y=meanAbsDifferenceSampleVSgoldstandard, group=1)) +
  geom_point() +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  labs(x="Input DNA",y="Difference between Sample and Gold Standard",
       title="QC - Differences between Sample and Gold Standard") 

ggplot(data=summary, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-23,23) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age") 

young_subset <- subset_data[subset_data$age <= 50,]
old_subset <- subset_data[subset_data$age > 50,]

summary_young <- getSummary(young_subset,"diff", "type")
summary_young$type <-  c("Naive","Central Memory","Effector Memory","TEMRA")
summary_young$type <-factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=summary_young, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-23,23) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age - < 50 years old") 

summary_old <- getSummary(old_subset,"diff", "type")
summary_old$type <- c("Naive","Central Memory","Effector Memory","TEMRA")
summary_old$type <- factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=summary_old, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-23,23) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age - >50 years old") 

ggplot(data=subset_data,aes(x=age, y=DNAmAge, color=type)) +
  geom_point() +
  theme_classic() +
  xlim(10,80) + ylim(10,80) +
  theme(text = element_text(size = 15)) +
  geom_abline(linetype="dotted") +
  geom_hline(yintercept = 0,linetype="dotted") +
  labs(x="Donor Age",y="Predicted Age", title="Age vs. Predicted Age Per Donor") 

#percent_meth <- readRDS("./ResultsClock/rgset.RDS")
#colData(percent_meth)$Sample_Well
#MsetEx <- preprocessIllumina(percent_meth)
#GMsetEx <- mapToGenome(MsetEx)
#qc <- minfiQC(GMsetEx, fixOutliers = TRUE, verbose = TRUE)
#GMsetEx <- qc$object
#ratio <- ratioConvert(GMsetEx)
#beta_values <- getBeta(ratio)
#beta_values <- data.frame(na.omit(beta_values))
#colnames(beta_values) <- colData(percent_meth)$Sample_Well
#beta_values<-beta_values[,all_data$SampleID]

#plotQC(qc$qc)

#all_data$old <- all_data$age >65

#annotations <- read.delim("~/Desktop/Data/subset_data/EPIC.hg38.manifest.tsv")
#matched_positions <- match(rownames(beta_values),annotations$probeID) #Finds valid matches in table.
#matched_symbols <- annotations$gene[matched_positions]
#cpg_data <- beta_values
#write.csv(cpg_data,"cpg_data.csv")
#beta_values$gene <- matched_symbols
#beta_values <- beta_values[!is.na(beta_values$gene),]
#cpg_data <- cpg_data[!is.na(beta_values$gene),]
#mapping <- cbind(cpg_data,beta_values)
#mapping$row_name <- make.names(beta_values$gene,unique=TRUE)
#rownames(beta_values) <- make.names(beta_values$gene,unique=TRUE)
#beta_values <- beta_values[,1:31]
#beginning_positions <- match(mapping$cpg,annotations$probeID)
#mapping$begin <- annotations$CpG_beg[beginning_positions]
#mapping$end <- annotations$CpG_end[beginning_positions]
#ditch the sabgal samples
#beta_values <- beta_values[all_data$sabgal_sample==FALSE,]
#write.csv(mapping,"mapping.csv")
mapping <- read.csv("mapping.csv")
beta_values <- read.csv("~/Desktop/Data/subset_data/beta_values.csv", row.names=1)
beta_values <- beta_values[,colnames(beta_values) %in% all_data$SampleID]

#limma differential methylation analysis
celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45553","R45804","R45504",
                                              "R45741","R45805"))
sabgal_sample <- factor(all_data$sabgal_sample,levels=c("TRUE","FALSE"))
old_group <- factor(all_data$old,levels=c(TRUE,FALSE))
design <- model.matrix(~0 + donor_group + celltype_group)

fit.reduced <- lmFit(beta_values,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
diff_exp <-topTable(fit.reduced,coef=10,number=1000000)
diff_exp[order(diff_exp$logFC),]
diff_exp[order(diff_exp$logFC, decreasing=TRUE),]

diff_exp$color=factor(case_when(diff_exp$adj.P.Val < .05 & abs(diff_exp$logFC) >= .6 ~ "purple",
                                           (diff_exp$adj.P.Val < .05 & abs(diff_exp$logFC) < .6) ~ "red",
                                           (diff_exp$adj.P.Val > .05 & abs(diff_exp$logFC) >= .6) ~ "blue",
                                           (diff_exp$adj.P.Val > .05 & abs(diff_exp$logFC) < .6) ~ "gray"))
diff_exp$delabel <- NA
diff_exp$delabel[diff_exp$color=="purple"] <- rownames(diff_exp)[diff_exp$color=="purple"]
ggplot(data=diff_exp, aes(x=logFC, y=-log10(adj.P.Val), color=color)) + 
  geom_point() +
  xlim(-1,1) + ylim(0,20) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  scale_colour_identity() +
  labs(title="Volcano Plot - TEMRA vs. Naive")

beta_rotated <- t(beta_values)
umap <- umap(beta_rotated)
umap_plot_df <- data.frame(umap$layout) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::inner_join(all_data, by = "SampleID")
umap_plot_df$type <- as.character(umap_plot_df$type)
umap_plot_df$type <- factor(umap_plot_df$type, levels=c("naive", "central_memory", "effector_memory","temra"))
ggplot(
  umap_plot_df,
  aes(x = X1, y = X2, color=type)) +
  labs(x="UMAP Component 1", y="UMAP Component 2", title = "UMAP Visualization") +
  theme_classic() +
  geom_point() # Plot individual points to make a scatterplot

ggplot(umap_plot_df,aes(x=type,y=-X1, color=age)) + 
  theme_classic() +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#messing with pseudotime

Y <- log2(beta_values+ 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
Y <- Y[var1K, ]  # only counts for variable genes

t <- -umap_plot_df$X1
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:150]  
heatdata <- as.matrix(beta_values[rownames(beta_values) %in% topgenes, order(t, na.last = NA)])
heatclus <- umap_plot_df$type[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', cexRow = 1.5, fontsize = 15)

write.csv(gsub("\\..*","",topgenes),"list.csv")

Y <- log2(beta_values+ 1)
var25K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:25000]
Y <- Y[var25K, ]  # only counts for variable genes
#Let's look at HDAC4
HDAC4 <- Y[rownames(Y)[str_detect(rownames(Y), "HDAC4")],]
HDAC4 <- colSums(HDAC4)/nrow(data.frame(HDAC4))
all_data$HDAC4 <- HDAC4
ggplot(all_data,aes(x=type,y=HDAC4)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at IFNG
IFNG <- Y[rownames(Y)[str_detect(rownames(Y), "IFNG")],]
IFNG <- colSums(IFNG)/nrow(data.frame(IFNG))
all_data$IFNG <- IFNG
ggplot(all_data,aes(x=type,y=IFNG)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at CCR2
CCR2 <- Y[rownames(Y)[str_detect(rownames(Y), "CCR2")],]
CCR2 <- colSums(CCR2)/nrow(data.frame(CCR2))
all_data$CCR2 <- CCR2
ggplot(all_data,aes(x=type,y=CCR2)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at ZBTB38
ZBTB38 <- Y[rownames(Y)[str_detect(rownames(Y), "ZBTB38")],]
ZBTB38 <- colSums(ZBTB38)/nrow(data.frame(ZBTB38))
all_data$ZBTB38 <- ZBTB38
ggplot(all_data,aes(x=type,y=ZBTB38)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at FLT3
FLT3 <- Y[rownames(Y)[str_detect(rownames(Y), "FLT3")],]
FLT3 <- colSums(FLT3)/nrow(data.frame(FLT3))
all_data$FLT3 <- FLT3
ggplot(all_data,aes(x=type,y=FLT3)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at TOX3
TOX3 <- Y[rownames(Y)[str_detect(rownames(Y), "TOX3")],]
TOX3 <- colSums(TOX3)/nrow(data.frame(TOX3))
all_data$TOX3 <- TOX3
ggplot(all_data,aes(x=type,y=TOX3)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="% Methylated",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at FAS
FAS <- Y[rownames(Y)[str_detect(rownames(Y), "FAS")],]
FAS <- colSums(FAS)/nrow(data.frame(FAS))
all_data$FAS <- FAS
ggplot(all_data,aes(x=type,y=FAS)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="% Methylated",title="FAS decreases in DNA methylation as differentiation proceeds") +
  geom_point()

#Let's look at CD28
CD28 <- Y[rownames(Y)[str_detect(rownames(Y), "CD28")],]
CD28 <- colSums(CD28)/nrow(data.frame(CD28))
all_data$CD28 <- CD28
ggplot(all_data,aes(x=type,y=CD28)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="% Methylated",title="CD28 increases in DNA methylation as differentiation proceeds") +
  geom_point()

#Let's look at LAG3
LAG3 <- Y[rownames(Y)[str_detect(rownames(Y), "LAG3")],]
LAG3 <- colSums(LAG3)/nrow(data.frame(LAG3))
all_data$LAG3 <- LAG3
ggplot(all_data,aes(x=type,y=LAG3)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at HAVCR2
HAVCR2 <- Y[rownames(Y)[str_detect(rownames(Y), "HAVCR2")],]
HAVCR2 <- colSums(HAVCR2)/nrow(data.frame(HAVCR2))
all_data$HAVCR2 <- HAVCR2
ggplot(all_data,aes(x=type,y=HAVCR2)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at TNFRSF9
TNFRSF9 <- Y[rownames(Y)[str_detect(rownames(Y), "TNFRSF9")],]
TNFRSF9 <- colSums(TNFRSF9)/nrow(data.frame(TNFRSF9))
all_data$TNFRSF9 <- TNFRSF9
ggplot(all_data,aes(x=type,y=TNFRSF9)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at SORL1
SORL1 <- Y[rownames(Y)[str_detect(rownames(Y), "SORL1")],]
SORL1 <- colSums(SORL1)/nrow(data.frame(SORL1))
all_data$SORL1 <- SORL1
ggplot(all_data,aes(x=type,y=SORL1)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at DNMT3A
DNMT3A <- Y[rownames(Y)[str_detect(rownames(Y), "DNMT3A")],]
DNMT3A <- colSums(DNMT3A)/nrow(data.frame(DNMT3A))
all_data$DNMT3A <- DNMT3A
ggplot(all_data,aes(x=type,y=DNMT3A)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at RASA2
RASA2 <- Y[rownames(Y)[str_detect(rownames(Y), "RASA2")],]
RASA2 <- colSums(RASA2)/nrow(data.frame(RASA2))
all_data$RASA2 <- RASA2
ggplot(all_data,aes(x=type,y=RASA2)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at GZMB
GZMB <- Y[rownames(Y)[str_detect(rownames(Y), "GZMB")],]
GZMB <- colSums(GZMB)/nrow(data.frame(GZMB))
all_data$GZMB <- GZMB
ggplot(all_data,aes(x=type,y=GZMB)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()
  
beta_values["IL7R",]
Y["HAVCR2",]
Y_down <- Y[(Y$A1 < Y$A3) & (rownames(Y) %in% topgenes),]
write.csv(gsub("\\..*","",rownames(Y_down)),"more_exp_in_temra.csv")

Y_up <- Y[(Y$A1 > Y$A3) & (rownames(Y) %in% topgenes),]
write.csv(gsub("\\..*","",rownames(Y_up)),"more_exp_in_naive.csv")

#Let's focus on HDAC4 specifically.
all_hdac4 <- diff_exp[rownames(diff_exp)[str_detect(rownames(diff_exp), "HDAC4")],]
all_hdac4$site <- rownames(all_hdac4)
match_starts <- match(rownames(all_hdac4),mapping$row_name)
all_hdac4$start <-mapping$begin[match_starts]-239233056
all_hdac4$position <- case_when(all_hdac4$start < 0 ~ "upstream",
                                all_hdac4$start > 0 & all_hdac4$start < 10000 ~ "in gene",
                                all_hdac4$start > 10000 ~ "downstream")
ggplot(all_hdac4,aes(x=(start),y=logFC,color=position) ) +
  theme_classic() + geom_bar(stat="identity",position=position_dodge(width=50),size=1.2) +
  xlab("Position Relative to HDAC4 ATG") +
  ylim(-1,1) + xlim(-7500,4000) +
  geom_vline(xintercept=0,linetype="dotted") + geom_vline(xintercept=9000,linetype="dotted") +
  geom_hline(yintercept=0) +
  ylab("DNAm logFC") +
  theme(text = element_text(size = 12))    +
  ggtitle("HDAC4 Upstream Region Shows Differential Methylation Patterns between TEMRA and Naive Cells")

#Let's focus on DNMT3A specifically.
all_DNMT3A <- diff_exp[rownames(diff_exp)[str_detect(rownames(diff_exp), "DNMT3A")],]
all_DNMT3A$site <- rownames(all_DNMT3A)
match_starts <- match(rownames(all_DNMT3A),mapping$row_name)
all_DNMT3A$start <-mapping$begin[match_starts]-25342590
all_DNMT3A$position <- case_when(all_DNMT3A$start < 0 ~ "upstream",
                                 all_DNMT3A$start > 0 & all_DNMT3A$start < 10000 ~ "in gene",
                                 all_DNMT3A$start > 10000 ~ "downstream")
ggplot(all_DNMT3A,aes(x=(start),y=logFC,color=position) ) +
  theme_classic() + geom_bar(stat="identity",position=position_dodge(width=50),size=1.2) +
  xlab("Position Relative to DNMT3A ATG") +
  ylim(-1,1) + xlim(-7500,4000) +
  geom_vline(xintercept=0,linetype="dotted") +
  geom_hline(yintercept=0) +
  ylab("DNAm logFC") +
  theme(text = element_text(size = 12))    +
  ggtitle("DNMT3A Upstream Region Shows Differential Methylation Patterns between TEMRA and Naive Cells")

#Let's focus on TET2 specifically.
CD28 <- diff_exp[rownames(diff_exp)[str_detect(rownames(diff_exp), "CD28")],]
CD28$site <- rownames(CD28)
match_starts <- match(rownames(CD28),mapping$row_name)
CD28$start <-mapping$begin[match_starts]-203706509
CD28$position <- case_when(CD28$start < 0 ~ "upstream",
                           CD28$start > 0 & CD28$start < 10000 ~ "in gene",
                           CD28$start > 10000 ~ "downstream")
ggplot(CD28,aes(x=(start),y=logFC,color=position) ) +
  theme_classic() + geom_bar(stat="identity",position=position_dodge(width=50),size=1.2) +
  xlab("Position Relative to CD28 ATG") +
  ylim(-1,1) + xlim(-75000,4000) +
  geom_vline(xintercept=0,linetype="dotted") +
  geom_hline(yintercept=0) +
  ylab("DNAm logFC") +
  theme(text = element_text(size = 12))    +
  ggtitle("CD28 Upstream Region Shows Differential Methylation Patterns between TEMRA and Naive Cells")
