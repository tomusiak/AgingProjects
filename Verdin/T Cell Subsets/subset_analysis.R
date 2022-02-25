#Grabs some useful scripts.
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

#Sets location of data 
setwd("/home/atom/Desktop/Data/subset_data") #Sets directory.

#Libraries to import.
library("cgageR")
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
library(rstatix)
library(ggpubr)

#Read in clock data and metadata. Merge them and process them.
clock_data <- read.csv("~/Desktop/Data/subset_data/clock_data.csv")
subset_metadata <- read.csv("~/Desktop/Data/subset_data/subset_metadata.csv")
all_data <- merge(clock_data,subset_metadata)
all_data$type <- as.character(all_data$type)
all_data$type <- factor(all_data$type, levels=c("naive", "central_memory", "effector_memory","temra"))

# #Read in CpG data. Process and acquire QC. Move around columns to match metadata.
# percent_meth <- readRDS("./ResultsClock/rgset.RDS")
# MsetEx <- preprocessNoob(percent_meth)
# GMsetEx <- mapToGenome(MsetEx)
# qc <- minfiQC(GMsetEx, fixOutliers = TRUE, verbose = TRUE)
# GMsetEx <- qc$object
# ratio <- ratioConvert(GMsetEx)
# beta_values <- getBeta(ratio)
# m_values <- getM(ratio)
# beta_values <- data.frame(na.omit(beta_values))
# m_values <- data.frame(na.omit(m_values))
# colnames(beta_values) <- colData(percent_meth)$Sample_Well
# colnames(m_values) <- colData(percent_meth)$Sample_Well
# beta_values<-beta_values[,all_data$SampleID]
# m_values<-m_values[,all_data$SampleID]
# 
# #Read in annotations to create a 'mapping table' that links together metadata, CpG data, and 
# #clock information.
# annotations <- read.delim("~/Desktop/Data/subset_data/EPIC.hg38.manifest.tsv")
# matched_positions <- match(rownames(beta_values),annotations$probeID) #Finds valid matches in table.
# matched_symbols <- annotations$gene[matched_positions]
# beta_values$gene <- matched_symbols
# m_values$gene <- matched_symbols
# beta_values <- na.omit(beta_values)
# m_values <- na.omit(m_values)
# cpg_data <- beta_values
# write.csv(cpg_data,"cpg_data.csv")
# unique_names <- make.names(beta_values$gene,unique=TRUE)
# rownames(beta_values) <- unique_names
# rownames(m_values) <- unique_names
# mapping <- cbind(cpg_data,beta_values)
# mapping$row_name <-rownames(beta_values) 
# beta_values <- beta_values[,1:32]
# m_values <- m_values[1:32]
# mapping$cpg <- rownames(mapping)
# beginning_positions <- match(mapping$cpg,annotations$probeID)
# mapping$begin <- annotations$CpG_beg[beginning_positions]
# mapping$end <- annotations$CpG_end[beginning_positions]
# write.csv(beta_values,"beta_values.csv")
# write.csv(mapping,"mapping.csv")
# write.csv(m_values,"m_values.csv")
beta_values <- read.csv("~/Desktop/Data/subset_data/beta_values.csv", row.names=1)
mapping <- read.csv("~/Desktop/Data/subset_data/mapping.csv", row.names=1)
m_values <- read.csv("~/Desktop/Data/subset_data/m_values.csv", row.names=1)

#Filter out unwanted data from metadata, mapping, and beta values.
keep <- all_data$SampleID[all_data$donor != "R45553" & all_data$sabgal_sample == FALSE]
all_data <- all_data[all_data$SampleID %in% keep,]
beta_values <- beta_values[,colnames(beta_values) %in% keep]
beta_values <- beta_values[order(all_data$type)]
m_values <- m_values[,colnames(m_values) %in% keep]
m_values <- m_values[order(all_data$type)]
all_data <- all_data[order(all_data$type),]

#Calculate epigenetic clock acceleration
all_data$diff <- all_data$DNAmAge-all_data$age
summary <- getSummary(all_data,"diff", "type")
summary$type <- c("Naive","Central Memory","Effector Memory","TEMRA")
summary$type <- factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

#subset ages
naive_data <- all_data[all_data$type == "naive",]
cm_data <- all_data[all_data$type == "central_memory",]
em_data <- all_data[all_data$type == "effector_memory",]
temra_data <- all_data[all_data$type == "temra",]

#t tests
t.test(em_data$IEAA,naive_data$IEAA,paired=TRUE)
t_test_data <- data.frame(all_data[all_data$donor != "R45553",])
t_test_data <- t_test_data[c("IEAA","type","donor")]
t_test_data$type <- as.factor(t_test_data$type)
t_test_data$donor <- as.factor(t_test_data$donor)
test <- t_test_data %>% pairwise_t_test(IEAA ~ type, paired=TRUE,correction="bonferroni")  %>%
                                                          select(-df, -statistic, -p)
res.aov <- t_test_data %>% anova_test(IEAA ~ type) 
res.aov

#QC check for DNA quantity
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
  ylim(-20,20) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age") 

summary_IEAA <- getSummary(all_data,"IEAA", "type")

ggplot(data=summary_IEAA, aes(x=type, y=IEAA, group=1)) +
  geom_point() +
  theme_classic() +
  ylim(-5,5) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=IEAA-se, ymax=IEAA+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="IEAA", title="CD8+ T Cell Subset Differences With IEAA") 

summary_IEAA_donor <- getSummary(all_data,"IEAA", "donor")

ggplot(data=all_data, aes(x=donor, y=IEAA, group=1,color=type)) +
  geom_point() +
  theme_classic() +
  ylim(-10,10) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  labs(x="Donor",y="IEAA", title="CD8+ IEAA Differences Between Donors") 

young_subset <- all_data[all_data$age <= 50,]
old_subset <- all_data[all_data$age > 50,]

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

ggplot(data=all_data,aes(x=age, y=DNAmAge, color=type)) +
  geom_point() +
  theme_classic() +
  xlim(10,80) + ylim(10,80) +
  theme(text = element_text(size = 15)) +
  geom_abline(linetype="dotted") +
  geom_hline(yintercept = 0,linetype="dotted") +
  labs(x="Donor Age",y="Predicted Age", title="Age vs. Predicted Age Per Donor") 

#limma differential methylation analysis
celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45804","R45504",
                                              "R45741","R45805"))
sabgal_sample <- factor(all_data$sabgal_sample,levels=c("TRUE","FALSE"))
old_group <- factor(all_data$old,levels=c(TRUE,FALSE))
design <- model.matrix(~0 + donor_group + celltype_group)

fit.reduced <- lmFit(beta_values,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
diff_exp <-topTable(fit.reduced,coef=10,number=10000)

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

ggplot(umap_plot_df,aes(x=type,y=X2, color=age)) + 
  theme_classic() +
  labs(x="CD8 Cell Subtype", y="UMAP Component 2",title="UMAP Component 2 Tracks Cell Lineage") +
  geom_point()

#messing with pseudotime

Y <- beta_values
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
Y <- Y[var1K, ]  # only counts for variable genes
t <- umap_plot_df$X2
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:150]  
heatdata <- as.matrix(beta_values[rownames(beta_values) %in% topgenes, order(t, na.last = NA)])
heatclus <- umap_plot_df$type[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", 
                               visualizeData = 'transformed', cexRow = 1.5, fontsize = 15)
write.csv(gsub("\\..*","",topgenes),"list.csv")

var25K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:25000]
Y <- Y[var25K, ]  # only counts for variable genes

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

data(HorvathLongCGlist)
clock_list <- HorvathLongCGlist
matching_cpgs <- match(clock_list$MR_var,mapping$cpg)
gene_names <- mapping$row_name[matching_cpgs]
clock_changes <- na.omit(diff_exp[gene_names,1:6])
clock_changes <- clock_changes[order(clock_changes$adj.P.Val),]
clock_changes$color=factor(case_when(clock_changes$adj.P.Val < .05 & abs(clock_changes$logFC) >= .25 ~ "purple",
                                (clock_changes$adj.P.Val < .05 & abs(clock_changes$logFC) < .25) ~ "red",
                                (clock_changes$adj.P.Val > .05 & abs(clock_changes$logFC) >= .25) ~ "blue",
                                (clock_changes$adj.P.Val > .05 & abs(clock_changes$logFC) < .25) ~ "gray"))
clock_changes$delabel <- NA
clock_changes$delabel[clock_changes$color=="purple"] <- rownames(clock_changes)[clock_changes$color=="purple"]
ggplot(data=clock_changes, aes(x=logFC, y=-log10(adj.P.Val), color=color,label=delabel)) + 
  geom_point() +
  xlim(-1,1) + ylim(0,12) +
  theme_classic(base_size=15)  +
  geom_vline(xintercept=.25,linetype="dotted") +
  geom_vline(xintercept=-.25,linetype="dotted") +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_text(nudge_y=.2) +
  scale_colour_identity() +
  labs(title="Volcano Plot of Clock CpG Sites Changing between Naive & TEMRA Samples")
head(clock_changes)

#Let's look at differential methylation by subtype by gene.

sasp_genes <- c("IL6R","IL4R","NFKB","IL1A")
naive_genes <- c("CD45RA","CCR7","CD62L","CD127")
effector_genes <- c("CD45RO","GRZMB","CD69","TNFR1")
exhausted_genes <- c("TIGIT","LAG3","HAVCR2","TOX")
