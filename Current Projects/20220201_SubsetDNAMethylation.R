#Grabs some useful scripts.
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

#Sets location of data 
setwd("/home/atom/Desktop/Data/subset_data") #Sets directory.
#Libraries to import.
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(methylGSA)
library(plyr)
library(reshape2)
library("cgageR")
library(WGCNA)
library(limma)
library(dplyr)
require("biomaRt")
library("RColorBrewer")
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
library(gplots)
library(minfi)

#Read in clock data and metadata. Merge them and process them.
clock_data <- read.csv("~/Desktop/Data/subset_data/clock_data.csv")
subset_metadata <- read.csv("~/Desktop/Data/subset_data/subset_metadata.csv")
all_data <- merge(clock_data,subset_metadata)
all_data$type <- as.character(all_data$type)
all_data$type <- factor(all_data$type, levels=c("naive", "central_memory", "effector_memory","temra"))

# #Read in CpG data. Process and acquire QC. Move around columns to match metadata.
#percent_meth <- readRDS("./ResultsClock/rgset.RDS")
#MsetEx <- preprocessNoob(percent_meth)
#GMsetEx <- mapToGenome(MsetEx)
#qc <- minfiQC(GMsetEx, fixOutliers = TRUE, verbose = TRUE)
#GMsetEx <- qc$object
#ratio <- ratioConvert(GMsetEx)
#beta_values <- getBeta(ratio)
#m_values <- getM(ratio)
#beta_values <- data.frame(na.omit(beta_values))
#m_values <- data.frame(na.omit(m_values))
#colnames(beta_values) <- colData(percent_meth)$Sample_Well
#colnames(m_values) <- colData(percent_meth)$Sample_Well
#beta_values<-beta_values[,all_data$SampleID]
#m_values<-m_values[,all_data$SampleID]
#write.csv(beta_values,"beta_values.csv")
#write.csv(m_values,"m_values.csv")
beta_values <- read.csv("~/Desktop/Data/subset_data/beta_values.csv", row.names=1)
#m_values <- read.csv("~/Desktop/Data/subset_data/m_values.csv", row.names=1)

#Filter out unwanted data from metadata, mapping, and beta values.
keep <- all_data$SampleID[all_data$SampleID != "D3" & all_data$sabgal_sample == FALSE]
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

beta_rotated <- data.frame(beta_rotated)
all_data$differentiation <- umap_plot_df$X2

#limma differential methylation analysis
celltype_group <- factor(all_data$type,levels=c("naive","central_memory","effector_memory","temra"))
donor_group <- factor(all_data$donor,levels=c("R45690","R45740","R45804","R45504",
                                              "R45553","R45741","R45805"))
sabgal_sample <- factor(all_data$sabgal_sample,levels=c("TRUE","FALSE"))
old_group <- factor(all_data$old,levels=c(TRUE,FALSE))
diff_group <- all_data$differentiation
design <- model.matrix(~0 + donor_group + diff_group)

fit.reduced <- lmFit(beta_values,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
diff_exp <-topTable(fit.reduced,coef=8,number=1000000)
#messing with pseudotime
diff_exp_order <- diff_exp[order(diff_exp$adj.P.Val),]
pvals <- as.vector(diff_exp_order[,5])
names(pvals) <- rownames(diff_exp_order)
adjustment <- methylglm(pvals, array.type = "EPIC")

Y <- drop_na(beta_values[rownames(diff_exp_order), ]) 
t <- umap_plot_df$X2
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

topgenes <- names(sort(gam.pval, decreasing = TRUE))[1:200]  
heatdata <- as.matrix(beta_values[rownames(beta_values) %in% topgenes,
                                  order(t, na.last = NA)])
heatclus <- umap_plot_df$type[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus)
clusterExperiment::plotHeatmap(ce, clusterSamplesData = "dendrogramValue",nFeatures=150, 
                               visualizeData = 'transformed', cexRow = 1.5, fontsize = 15)

#Let's create a custom CpG annotation table.
genome_annotation <- read.delim("~/Desktop/Data/gencode.v39.annotation.gtf", header=FALSE, comment.char="#")

#Let's re-do the analysis using only the mapped CpGs.

#variable_cpgs <- names(sort(apply(beta_values, 1, var),decreasing = TRUE))[1:400000]
#variable_betas <- beta_values[variable_cpgs,]
#variable_annotation <- cpg_annotation[cpg_annotation$probeID %in% variable_cpgs,]

#cpg_table <- createCPGTable(variable_annotation,genome_annotation,1000,2500)
#write.csv(cpg_table,"cpg_table.csv")

cpg_table <- read.csv("~/Desktop/Data/subset_data/cpg_table.csv", row.names=1)

#Let's look at differential methylation by subtype by gene
dna_modifiers <- c("DNMT1","DNMT3A","DNMT3B",
                   "TET1","TET2","TET3",
                   "APOBEC3A","APOBEC3D","APOBEC3G","APOBEC3H")

histone_modifiers <- c("HDAC1","HDAC4",
                       "HDAC6","HDAC9", "HDAC11",
                       "SIRT2","SIRT4","SIRT5",
                       "SIRT6","SIRT7", "KAT2A","KAT2B",
                       "KAT6A","KAT6B",
                       "EHMT1","EHMT2",
                       "NSD1","PRDM2","SET","SETD3","SETD5","SETD6",
                       "SETMAR","SMYD1","SMYD2","SMYD3","SMYD5",
                       "SUV39H1","SUV39H2","ATR")

interesting_histones <- c("SET","HDAC1","KAT2A",
                          "HDAC11","HDAC9","SIRT2","KAT6A",
                          "KAT6B","KAT2B")

interleukins <- c("IL3","IL4","IL5","IL6",
                  "IL9","IL10","IL11",
                  "IL13","IL15","IL16","IL17D","IL17A",
                  "IL17B","IL17C","IL25")

known_diff_TFs <- c("TCF7","EOMES", "FOXO1",
                    "BATF","IRF4",
                    "TBX21","RUNX3","PRDM1",
                    "LEF1")

transcription_factors <- as.list(read.table("~/Desktop/Data/subset_data/TFs.txt", 
                                            quote="\"", comment.char=""))$V1

dna <- getDiffMethylationList(dna_modifiers,cpg_table)
rownames(dna) <- c("Naive","CM","EM","TEMRA")
dna <- dna[,-1]
dna <- log2(dna / rbind(dna[1,],dna[1,],dna[1,],dna[1,]))
colors = colorRampPalette(c("blue", "black", "green"))(n = 30)

histone <- getDiffMethylationList(histone_modifiers,cpg_table)
rownames(histone) <- c("Naive","CM","EM","TEMRA")
histone <- histone[,-1]
histone <- log2(histone / rbind(histone[1,],histone[1,],histone[1,],histone[1,]))

ils <- getDiffMethylationList(interleukins,cpg_table)
rownames(ils) <- c("Naive","CM","EM","TEMRA")
ils <- ils[,-1]
ils <- log2(ils / rbind(ils[1,],ils[1,],ils[1,],ils[1,]))

TFs <- getDiffMethylationList(transcription_factors,cpg_table)
TFs <- TFs %>% select_if(~all(!is.na(.)))
rownames(TFs) <- c("Naive","CM","EM","TEMRA")
TFs <- TFs[,-1]
TFs <- log2(TFs / rbind(TFs[1,],TFs[1,],TFs[1,],TFs[1,]))

heatmap.2(t(dna),density.info="none",Rowv = TRUE,Colv=FALSE,dendrogram="row",trace="none",
          main="DNA Modifying Enzymes", revC=TRUE, 
          colsep=1:4, xlab="Cell Type", key.xlab="Diff Methylation (relative to Naive)",
          breaks=seq(-1.5,1.5,0.1),col=colors,symkey=F,
          margins =c(7,7),cexRow=1,cexCol=1.2)

heatmap.2(t(histone),density.info="none",Rowv = TRUE,dendrogram="row",trace="none",
          Colv=FALSE, 
          main="Histone Editing Enzymes", revC=TRUE,
          colsep=1:10, xlab="Cell Type", key.xlab="Diff Methylation (relative to Naive)", 
          breaks=seq(-1.5,1.5,0.1), col=colors,symkey=F, 
          margins =c(10,10),cexRow=.5,cexCol=1.2)

heatmap.2(t(ils),density.info="none",dendrogram="row", trace="none", Rowv=TRUE,
          Colv=FALSE, colsep=1:4,
          main="Interleukins",  revC=TRUE, 
          xlab="Cell Type", key.xlab="Diff Methylation (relative to Naive)", 
          breaks=seq(-1.5,1.5,0.1), col=colors, symkey=F,
          margins =c(10,10),cexRow=.8,cexCol=1.2)

heatmap.2(t(TFs),density.info="none",Rowv = TRUE,dendrogram="row", trace="none",
          Colv=FALSE, colsep=1:4,
          main="Transcription Factors",  revC=TRUE, 
          xlab="Cell Type", key.xlab="Diff Methylation (relative to Naive)", 
          breaks=seq(-3,3,0.2), col=colors, symkey=F,
          margins =c(10,10),cexRow=.1,cexCol=1.2)

factors <- data.frame(t(TFs))

hyper <- head(factors[order(factors$EM,decreasing=TRUE),],12)
hypo <- tail(factors[order(factors$EM,decreasing=TRUE),],12)
cool_TFs <- data.frame(t(rbind(hyper,hypo)))

heatmap.2(t(cool_TFs),density.info="none",Rowv = TRUE,dendrogram="row", trace="none",
          Colv=FALSE, colsep=1:4,
          main="Transcription Factors",  revC=TRUE, 
          xlab="Cell Type", key.xlab="Diff Methylation (relative to Naive)", 
          breaks=seq(-3,3,0.2), col=colors, symkey=F,
          margins =c(10,10),cexRow=.6,cexCol=1.2)

hyper_list<-getDiffMethylationList2(rownames(hyper),cpg_table)
melted_hyper<-drop_na(melt(hyper_list))
melted_hyper$type <- factor(melted_hyper$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(melted_hyper,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Cell Type") +
  ylab("Methylated %") +
  ggtitle("Most-Hypermethylated TFs with Differentiation") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

hypo_list<-getDiffMethylationList2(rownames(hypo),cpg_table)
melted_hypo<-drop_na(melt(hypo_list))
melted_hypo$type <- factor(melted_hypo$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(melted_hypo,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") + 
  xlab("Cell Type") +
  ylab("Methylated %") +
  ggtitle("Most-Hypomethylated TFs with Differentiation") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

dna_binder_list<-getDiffMethylationList2(dna_modifiers[1:9],cpg_table)
melted_dna<-drop_na(melt(dna_binder_list))
melted_dna$type <- factor(melted_dna$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(melted_dna,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  ggtitle("DNA-Modifying Enzymes") +
  xlab("Cell Type") +
  ylab("Methylated %") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

histones_list<-getDiffMethylationList2(interesting_histones,cpg_table)
melted_histones<-drop_na(melt(histones_list))
melted_histones$type <- factor(melted_histones$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(melted_histones,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Cell Type") +
  ylab("Methylated %") +
  ggtitle("Histone-Modifying Enzymes") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

complete_table <- getAllGeneMethylation(cpg_table)
complete_table <- na.omit(complete_table)
write.csv(complete_table,"complete_table.csv")

complete_table <- read.csv("~/Desktop/Data/subset_data/complete_table.csv", row.names=1)

#Let's look at UMAP per genes.
gene_table_rotated <- t(complete_table)
gene_umap <- umap(gene_table_rotated)
gene_umap_plot_df <- data.frame(gene_umap$layout) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::inner_join(all_data, by = "SampleID")
gene_umap_plot_df$type <- as.character(gene_umap_plot_df$type)
gene_umap_plot_df$type <- factor(gene_umap_plot_df$type, levels=c("naive", "central_memory", "effector_memory","temra"))
ggplot(
  gene_umap_plot_df,
  aes(x = X1, y = X2, color=type)) +
  labs(x="UMAP Component 1", y="UMAP Component 2", title = "UMAP Visualization - Genes") +
  theme_classic() +
  geom_point() # Plot individual points to make a scatterplot

ggplot(gene_umap_plot_df,aes(x=type,y=X2, color=age)) + 
  theme_classic() +
  labs(x="CD8 Cell Subtype", y="UMAP Component 2",title="UMAP Component 2 Tracks Cell Lineage") +
  geom_point()

all_data$differentiation <- -gene_umap_plot_df$X2
diff_group <- all_data$differentiation
design <- model.matrix(~0 + donor_group + diff_group)

fit.reduced_2 <- lmFit(complete_table,design)
fit.reduced_2 <- eBayes(fit.reduced_2, robust=TRUE)
summary(decideTests(fit.reduced_2))
diff_exp_2 <-topTable(fit.reduced_2,coef=8,number=400000)
diff_exp_order_2 <- diff_exp_2[order(diff_exp_2$adj.P.Val),]

diff_exp_order_2 <-
  diff_exp_order_2[!grepl("ENSG",rownames(diff_exp_order_2)),]

diff_exp_order_2 <-
  diff_exp_order_2[!grepl("MIR",rownames(diff_exp_order_2)),]

diff_exp_order_2 <-
  diff_exp_order_2[!grepl("RNU",rownames(diff_exp_order_2)),]

diff_exp_order_2 <-
  diff_exp_order_2[!grepl("LINC",rownames(diff_exp_order_2)),]

top <- rownames(head(diff_exp_order_2,12))

top_list<-getDiffMethylationList2(top,cpg_table)
top_list<-drop_na(melt(top_list))
top_list$type <- factor(top_list$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(top_list,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Cell Type") +
  ylab("Methylated %") +
  ggtitle("Most Differentially-Methylated Genes With Differentiation") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

#Creates a volcano plot using the corrected counts.
diff_exp_order_2$color=factor(case_when(diff_exp_order_2$adj.P.Val < .05 & abs(diff_exp_order_2$logFC) >= .10 ~ "purple",
                                        (diff_exp_order_2$adj.P.Val < .05 & abs(diff_exp_order_2$logFC) < .10) ~ "red",
                                        (diff_exp_order_2$adj.P.Val >= .05 & abs(diff_exp_order_2$logFC) >= .10) ~ "blue",
                                        (diff_exp_order_2$adj.P.Val >= .05 & abs(diff_exp_order_2$logFC) < .10) ~ "gray"))
diff_exp_order_2$delabel <- NA
diff_exp_order_2$delabel[diff_exp_order_2$color=="purple"] <- rownames(diff_exp_order_2)[diff_exp_order_2$color=="purple"]
ggplot(data=diff_exp_order_2, aes(x=logFC, y=-log10(adj.P.Val), color=color)) + 
  geom_point() +
  xlim(-.5,.5) + ylim(0,20) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .1,linetype="dotted") +
  geom_vline(xintercept = -.1,linetype="dotted") +
  scale_colour_identity() +
  labs(title="Volcano Plot for Differentiation")

#Let's look at aging!
all_data$older <- c("Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Younger","Older","Older","Older","Younger","Older",
                    "Younger","Older","Older","Older","Younger","Older")
age_group <- all_data$older
age_design <- model.matrix(~0 + celltype_group + age_group)

fit.reduced_age <- lmFit(complete_table,age_design)
fit.reduced_age <- eBayes(fit.reduced_age, robust=TRUE)
summary(decideTests(fit.reduced_age))
diff_exp_age <-topTable(fit.reduced_age,coef=5,number=400000)
diff_exp_order_age <- diff_exp_age[order(diff_exp_age$adj.P.Val),]

diff_exp_order_age <-
  diff_exp_order_age[!grepl("ENSG",rownames(diff_exp_order_age)),]
top_age <- rownames(head(diff_exp_order_age,12))

top_list_age<-getDiffMethylationListAge(top_age,cpg_table)
top_list_age<-drop_na(melt(top_list_age))

ggplot(top_list_age,aes(x=older,y=value,fill=older)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Age Group") +
  ylab("Methylated %") +
  ggtitle("Most Differentially-Methylated Genes - Aging") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

#Creates a volcano plot using the corrected counts.
diff_exp_order_age$color=factor(case_when(diff_exp_order_age$adj.P.Val < .05 & abs(diff_exp_order_age$logFC) >= .1 ~ "purple",
                                          (diff_exp_order_age$adj.P.Val < .05 & abs(diff_exp_order_age$logFC) < .1) ~ "red",
                                          (diff_exp_order_age$adj.P.Val >= .05 & abs(diff_exp_order_age$logFC) >= .1) ~ "blue",
                                          (diff_exp_order_age$adj.P.Val >= .05 & abs(diff_exp_order_age$logFC) < .1) ~ "gray"))
diff_exp_order_age$delabel <- NA
diff_exp_order_age$delabel[diff_exp_order_age$color=="purple"] <- rownames(diff_exp_order_age)[diff_exp_order_age$color=="purple"]
ggplot(data=diff_exp_order_age, aes(x=logFC, y=-log10(adj.P.Val), color=color)) + 
  geom_point() +
  xlim(-.5,.5) + ylim(0,20) +
  theme_classic(base_size=18)  +
  geom_hline(yintercept = 1.2,linetype="dotted") +
  geom_vline(xintercept = .1,linetype="dotted") +
  geom_vline(xintercept = -.1,linetype="dotted") +
  scale_colour_identity() +
  labs(title="Volcano Plot for Aging")

fit_aging <- lmFit(beta_values,age_design)
fit_aging <- eBayes(fit_aging, robust=TRUE)
summary(decideTests(fit_aging))
diff_exp_aging <-topTable(fit_aging,coef=5,number=100000)
#messing with pseudotime
diff_exp_aging_order <- diff_exp_aging[order(diff_exp_aging$adj.P.Val),]
pvals_aging <- as.vector(diff_exp_aging_order[,5])
names(pvals_aging) <- rownames(diff_exp_aging_order)
diff_aging <- methylglm(pvals_aging, array.type = "EPIC")

all_data$PC1 <- -gene_umap_plot_df$X1
PC1_group <- all_data$PC1
PC1_design <- model.matrix(~0 + donor_group + PC1_group)

fit_reduced_PC1 <- lmFit(complete_table,PC1_design)
fit_reduced_PC1 <- eBayes(fit_reduced_PC1, robust=TRUE)
summary(decideTests(fit_reduced_PC1))
diff_exp_PC1 <-topTable(fit_reduced_PC1,coef=8,number=400000)
diff_exp_order_PC1 <- diff_exp_PC1[order(diff_exp_PC1$adj.P.Val),]

diff_exp_order_PC1 <-
  diff_exp_order_PC1[!grepl("ENSG",rownames(diff_exp_order_PC1)),]

top_PC1 <- rownames(head(diff_exp_order_PC1,12))
top_list_PC1<-getDiffMethylationList2(top_PC1,cpg_table)
top_list_PC1$type <- factor(top_list_PC1$type,levels=c("Naive","CM","EM","TEMRA"))
top_list_PC1<-drop_na(melt(top_list_PC1))
ggplot(top_list_PC1,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("> 60 Years Old") +
  ylab("Methylated %") +
  ggtitle("Most Differentially-Methylated Genes with UMAP Component 1") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

#TFs going UP
TFs_hypermeth <- diff_exp_2[rownames(factors),]
TFs_hypermeth <- TFs_hypermeth[TFs_hypermeth$adj.P.Val < .01 & TFs_hypermeth$logFC < -.05,]
TFs_hypermeth <- TFs_hypermeth[order(TFs_hypermeth$logFC),]
top_TFs_hypermeth <- TFs_hypermeth[1:12,]

top_TFs_hypermeth<-getDiffMethylationList2(rownames(top_TFs_hypermeth),cpg_table)
top_TFs_hypermeth<-drop_na(melt(top_TFs_hypermeth))
top_TFs_hypermeth$type <- factor(top_TFs_hypermeth$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(top_TFs_hypermeth,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Cell Type") +
  ylab("Methylated %") +
  ggtitle("Most Hypomethylated TFs With Differentiation") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")


#TFs going DOWN
TFs_hypometh <- diff_exp_2[rownames(factors),]
TFs_hypometh <- TFs_hypometh[TFs_hypometh$adj.P.Val < .01 & TFs_hypometh$logFC > -.05,]
TFs_hypometh <- TFs_hypometh[order(TFs_hypometh$logFC,decreasing=TRUE),]
top_TFs_hypometh <- TFs_hypometh[1:12,]

top_TFs_hypometh<-getDiffMethylationList2(rownames(top_TFs_hypometh),cpg_table)
top_TFs_hypometh<-drop_na(melt(top_TFs_hypometh))
top_TFs_hypometh$type <- factor(top_TFs_hypometh$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(top_TFs_hypometh,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Cell Type") +
  ylab("Methylated %") +
  ggtitle("Most Hypermethylated TFs With Differentiation") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

#Canonical TFs
important_TFs <- diff_exp_2[known_diff_TFs,]
important_TFs<-getDiffMethylationList2(rownames(important_TFs),cpg_table)
important_TFs<-drop_na(melt(important_TFs))
important_TFs$type <- factor(important_TFs$type,levels=c("Naive","CM","EM","TEMRA"))

ggplot(important_TFs,aes(x=type,y=value,fill=type)) +
  facet_wrap( ~ name, nrow = 3) +
  geom_violin(alpha=.5) + theme_classic() +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5,
               position="dodge") +
  theme(legend.position = "none") +
  xlab("Cell Type") +
  ylab("Methylated %") +
  ggtitle("Most Hypermethylated TFs With Differentiation") +
  stat_summary(fun = "mean",
               geom = "pointrange",
               color = "violet") +
  scale_fill_brewer(palette="RdYlBu")

#Looking at clock CpGs
#data(HorvathLongCGlist)
#clock_list <- HorvathLongCGlist
#matching_cpgs <- match(clock_list$MR_var,mapping$cpg)
#gene_names <- mapping$row_name[matching_cpgs]
#clock_changes <- na.omit(diff_exp[gene_names,1:6])
#clock_changes <- clock_changes[order(clock_changes$adj.P.Val),]
#clock_changes$color=factor(case_when(clock_changes$adj.P.Val < .05 & abs(clock_changes$logFC) >= .3 ~ "purple",
#                                     (clock_changes$adj.P.Val < .05 & abs(clock_changes$logFC) < .3) ~ "red",
#                                     (clock_changes$adj.P.Val > .05 & abs(clock_changes$logFC) >= .3) ~ "blue",
#                                     (clock_changes$adj.P.Val > .05 & abs(clock_changes$logFC) < .3) ~ "gray"))
#clock_changes$delabel <- NA
#clock_changes$delabel[clock_changes$color=="purple"] <- rownames(clock_changes)[clock_changes$color=="purple"]
#ggplot(data=clock_changes, aes(x=logFC, y=-log10(adj.P.Val), color=color,label=delabel)) + 
#  geom_point() +
#  xlim(-1,1) + ylim(0,12) +
#  theme_classic(base_size=15)  + 
#  geom_vline(xintercept=.25,linetype="dotted") +
#  geom_vline(xintercept=-.25,linetype="dotted") +
#  geom_hline(yintercept = 1.2,linetype="dotted") +
#  geom_text(nudge_y=.2) +
#  scale_colour_identity() +
#  labs(title="Volcano Plot of Clock CpG Sites Changing between Naive & TEMRA Samples")
#head(clock_changes)

#let's test a change