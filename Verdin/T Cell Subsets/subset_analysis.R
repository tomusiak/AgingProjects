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
summary$type <- c("Central Memory","Effector Memory","Naive","TEMRA")
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

young_subset <- subset_data[subset_data$age <= 65,]
old_subset <- subset_data[subset_data$age > 65,]

summary_young <- getSummary(young_subset,"diff", "type")
summary_young$type <- c("Central Memory","Effector Memory","Naive","TEMRA")
summary_young$type <- factor(summary_young$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=summary_young, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-23,23) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age - < 65 years old") 

summary_old <- getSummary(old_subset,"diff", "type")
summary_old$type <- c("Central Memory","Effector Memory","Naive","TEMRA")
summary_old$type <- factor(summary_old$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=summary_old, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-23,23) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age - >65 years old") 


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
#beta_values$gene <- matched_symbols
#beta_values <- beta_values[!is.na(beta_values$gene),]
#rownames(beta_values) <- make.names(beta_values$gene,unique=TRUE)
#beta_values <- beta_values[,1:31]

#ditch the sabgal samples
#beta_values <- beta_values[all_data$sabgal_sample==FALSE,]

beta_values <- read.csv("~/Desktop/Data/subset_data/beta_values.csv", row.names=1)

corSample=cor(beta_values,use = "p")
hierSample=hclust(as.dist(1-corSample), method="a")
branch1=cutreeStatic(hierSample,.03,minSize=2)
all_data$ClusteringBranch=branch1
all_data$CellColor= labels2colors(all_data$type)
all_data$ClusteringColor=matchLabels(labels2colors(branch1), as.character(all_data$CellColor) )
datColors=data.frame(branch= all_data$ClusteringColor,
                     CellType=all_data$CellColor, 
                     Donor = labels2colors(all_data$donor),
                     Older65 = labels2colors(all_data$old)) 
plotDendroAndColors(hierSample, colors=datColors)

#limma differential methylation analysis
celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45553","R45804","R45504",
                                              "R45741","R45805"))
sabgal_sample <- factor(all_data$sabgal_sample,levels=c("TRUE","FALSE"))
old_group <- factor(all_data$old,levels=c(TRUE,FALSE))
design <- model.matrix(~celltype_group + old_group  + donor_group)

fit.reduced <- lmFit(beta_values,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
diff_exp <-topTable(fit.reduced,coef=3,number=1000000)
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

#let's imagine it all as a time course experiment
ages <- factor(all_data$age,c(30,60,34,70,73,29,77))
design <- model.matrix(~0+ages)
colnames(design) <- c("a30","a60","a34","a70","a73","a29","a77")
fit <- lmFit(beta_values, design)
cont.wt <- makeContrasts(
 "a73-a30","a60-a30",
 levels=design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")
 
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

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:200]  
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
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at FAS
FAS <- Y[rownames(Y)[str_detect(rownames(Y), "FAS")],]
FAS <- colSums(FAS)/nrow(data.frame(FAS))
all_data$FAS <- FAS
ggplot(all_data,aes(x=type,y=FAS)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
  geom_point()

#Let's look at CD28
CD28 <- Y[rownames(Y)[str_detect(rownames(Y), "CD28")],]
CD28 <- colSums(CD28)/nrow(data.frame(CD28))
all_data$CD28 <- CD28
ggplot(all_data,aes(x=type,y=CD28)) + 
  theme_classic() +
  ylim(0,1) +
  labs(x="CD8 Cell Subtype", y="UMAP Component 1",title="UMAP Component 1 Tracks Cell Lineage") +
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
