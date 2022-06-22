#Sets location of Olga's data.
source("AgingProjects/Useful Scripts/generally_useful.R")
setwd("Data/Olga") #Sets directory.

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
library("Mus.musculus")
library(minfi)
library(umap)
library(ggvenn)

#Loading in data.
load("./data/NormalizedData/all_probes_sesame_normalized.Rdata")
datSample=read.csv("./data/SampleSheetAgeN80version4.csv")
annotations <-  readRDS("./tools/geneAnnotations/AnnotationAminHaghani/Latest versions/Mus_musculus_wsbeij.WSB_EiJ_v1.100.Mingjia.AminV2.RDS")
dat0=data.frame(normalized_betas_sesame)
rownames(dat0) <- dat0$CGid
dat0 <- dat0[,2:81]
colnames(dat0) <- substring(colnames(dat0),2)
dat0 <- dat0[,datSample$Basename]
datSample$newID <- paste(datSample$ExternalSampleID,datSample$Genotype,sep="_")
colnames(dat0) <- datSample$newID
orig <- dat0

#Let's swap CG ID's for gene names, and remove any NA values.
matched_positions <- match(rownames(dat0),annotations$CGid) #Finds valid matches in table.
matched_symbols <- annotations$SYMBOL[matched_positions]
matched_category <- annotations$main_Categories[matched_positions]
dat0$Symbol <- matched_symbols
dat0$category <- matched_category
dat0 <- subset(dat0, dat0$category == "Exon" | dat0$category == "Promoter")
dat0 <- dat0[!is.na(dat0$Symbol),]
dat0 <- dat0[,1:(ncol(dat0)-2)]
# rownames(dat0) <- make.names(dat0$Symbol,unique=TRUE)
# dat0 <- dat0[,1:81]
# m_orig <- log2(orig/(1-orig))
# colnames(dat0)[-1]=paste0(datSample$OriginalOrderInBatch,datSample$Tissue)
# colnames(orig)=paste0(datSample$OriginalOrderInBatch,datSample$Tissue)

#Some QC checks and re-formatting performed by Horvath's group and copied here.
#QC primarily performed by cluster analysis
corSample=cor(dat0[,-ncol(dat0)], use="p")
hierSample=hclust(as.dist(1-corSample), method="a")
branch1=cutreeStatic(hierSample,.02,minSize=2)
datSample$ClusteringBranch=branch1
datSample$TissueColor= labels2colors(datSample$Tissue)
datSample$ClusteringColor=matchLabels(labels2colors(branch1), as.character(datSample$TissueColor) )
datColors=data.frame(branch= datSample$ClusteringColor,
                     Tissue=datSample$TissueColor, 
                     Genotype=labels2colors(datSample$Genotype),
                     CanBeUsedForAging=labels2colors(datSample$CanBeUsedForAgingStudies),
                     Female=ifelse(datSample$Female==1,"pink","lightblue") )
plotDendroAndColors(hierSample, colors=datColors)

#processing
data_no_symbols <- dat0[,-ncol(dat0)]
unique_names <- make.names(dat0$Symbol,unique=TRUE)
data_symbolized <- data_no_symbols
rownames(data_symbolized) <- unique_names
m_values <- log2(data_symbolized/(1-data_symbolized))

#let's combine CpG sites by gene..
aggregated <- aggregate(x = m_values,                # Specify data column
          by = list(dat0$Symbol),              # Specify group indicator
          FUN = mean)
rownames(aggregated) <- aggregated$Group.1
aggregated <- aggregated[,-1]

#PCA
rotated <- t(aggregated)
pca<-prcomp(rotated,scale=TRUE)
var_explained = pca$sdev^2 / sum(pca$sdev^2)
pca <- data.frame(pca$x)
pca$genotype <- datSample$Genotype
pca$tissue <- datSample$Tissue
pca$sex <- datSample$Sex
pca$age <- datSample$Age
ggplot(data=data.frame(pca),aes(x=PC1,y=PC8,color=genotype,shape=tissue)) +
  geom_point() + theme_classic()
qplot(c(1:20), var_explained[1:20]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1) + theme_classic()

#PCA for heart only.
rotated <- t(aggregated[grepl("heart",colnames(aggregated))])
pca<-prcomp(rotated,scale=TRUE)
var_explained = pca$sdev^2 / sum(pca$sdev^2)
pca <- data.frame(pca$x)
pca$genotype <- datSample[datSample$Tissue=="Heart",]$Genotype
pca$tissue <- datSample[datSample$Tissue=="Heart",]$Tissue
pca$sex <- datSample[datSample$Tissue=="Heart",]$Sex
pca$age <- datSample[datSample$Tissue=="Heart",]$Age
ggplot(data=data.frame(pca),aes(x=PC1,y=PC2,color=genotype,shape=sex)) +
  geom_point() + theme_classic()
qplot(c(1:10), var_explained[1:10]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1) + theme_classic()

#PCA for liver only.
rotated <- t(aggregated[grepl("liver",colnames(aggregated))])
pca<-prcomp(rotated,scale=TRUE)
var_explained = pca$sdev^2 / sum(pca$sdev^2)
pca <- data.frame(pca$x)
pca$genotype <- datSample[datSample$Tissue=="Liver",]$Genotype
pca$tissue <- datSample[datSample$Tissue=="Liver",]$Tissue
pca$sex <- datSample[datSample$Tissue=="Liver",]$Sex
pca$age <- datSample[datSample$Tissue=="Liver",]$Age
ggplot(data=data.frame(pca),aes(x=PC1,y=PC2,color=genotype,shape=sex)) +
  geom_point() + theme_classic()
qplot(c(1:10), var_explained[1:10]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1) + theme_classic()

#limma differential methylation analysis
genotype_group <- factor(datSample$Genotype,levels=c("WT","KO"))
tissue_group <- factor(datSample$Tissue,levels=c("Heart","Brain","Liver","Muscle","Cerebellum"))
sex_group <- factor(datSample$Sex)
design <- model.matrix(~ tissue_group + sex_group + genotype_group)
fit.reduced <- lmFit(aggregated,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
diff_exp <-topTable(fit.reduced,coef=7,number=30000)

#Differential variability analysis using missMethyl.
fitvar <- varFit(aggregated, design = design)
summary(decideTests(fitvar))
topDV <- topVar(fitvar, coef=7,number=300)

#Volcano plot
cols <- densCols(diff_exp$logFC, -log10(diff_exp$adj.P.Val),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= diff_exp$logFC, 
     y = -log10(diff_exp$adj.P.Val), 
     col=cols, panel.first=grid(),
     main="Volcano plot of KO vs. WT mice DNA methylation levels at a CpG site", 
     xlim=c(-1,1),
     ylim=c(0,12),
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     cex=0.6)
gn.selected <- abs(diff_exp$logFC) >.06 & (diff_exp$adj.P.Val) < .00001
text(diff_exp$logFC[gn.selected],
     -log10(diff_exp$adj.P.Val)[gn.selected],
     lab=rownames(diff_exp)[gn.selected ], cex=0.6)

#GO analysis
sig_dmrs <- data.frame( gsub("\\..*","",rownames(diff_exp[diff_exp$adj.P.Val <= .001 &
                                                            abs(diff_exp$logFC) >=.15,])))
not_sig_dmrs <- data.frame(gsub("\\..*","",rownames(diff_exp[diff_exp$adj.P.Val > .001&
                                                               abs(diff_exp$logFC) <.15,])))
colnames(sig_dmrs) <- "gene_ids"
colnames(not_sig_dmrs) <- "gene_ids"
sig_dmrs$is_candidate <- 1
not_sig_dmrs$is_candidate <- 0
combined <- rbind(sig_dmrs,not_sig_dmrs)
res_hyper_bg = go_enrich(combined, n_randsets=100,organism="Mus.musculus")
top_gos_hyper_bg = res_hyper_bg[[1]][1:5, 'node_id']
plot_stats = plot_anno_scores(res_hyper_bg, top_gos_hyper_bg)
plot_stats
res_hyper_bg[[1]][1:10, 'node_name']

#other stuff
sig_dmrs_full_data<-diff_exp[sig_dmrs$gene_ids,]
mean(diff_exp$logFC)

#let's subset - muscle only
muscle_only <- aggregated[grepl("muscle",colnames(aggregated))]
muscle_only_annot <- datSample[datSample$Tissue == "Muscle",]
genotype_group_muscle <- factor(muscle_only_annot$Genotype,levels=c("WT","KO"))
sex_group_muscle <- factor(muscle_only_annot$Sex)
design_muscle <- model.matrix(~ sex_group_muscle + genotype_group_muscle)
fit.reduced_muscle <- lmFit(muscle_only,design_muscle)
fit.reduced_muscle <- eBayes(fit.reduced_muscle, robust=TRUE)
summary(decideTests(fit.reduced_muscle))
diff_exp_muscle <-topTable(fit.reduced_muscle,coef=3,number=30000)

#let's subset - heart only
heart_only <- aggregated[grepl("heart",colnames(aggregated))]
heart_only_annot <- datSample[datSample$Tissue == "Heart",]
genotype_group_heart <- factor(heart_only_annot$Genotype,levels=c("WT","KO"))
sex_group_heart <- factor(heart_only_annot$Sex)
design_heart <- model.matrix(~ sex_group_heart + genotype_group_heart)
fit.reduced_heart <- lmFit(heart_only,design_heart)
fit.reduced_heart <- eBayes(fit.reduced_heart, robust=TRUE)
summary(decideTests(fit.reduced_heart))
diff_exp_heart <-topTable(fit.reduced_heart,coef=3,number=30000)

#let's subset - cerebellum only
cerebellum_only <- aggregated[grepl("cerebellum",colnames(aggregated))]
cerebellum_only_annot <- datSample[datSample$Tissue == "Cerebellum",]
genotype_group_cerebellum <- factor(cerebellum_only_annot$Genotype,levels=c("WT","KO"))
sex_group_cerebellum <- factor(cerebellum_only_annot$Sex)
design_cerebellum <- model.matrix(~ sex_group_cerebellum + genotype_group_cerebellum)
fit.reduced_cerebellum <- lmFit(cerebellum_only,design_cerebellum)
fit.reduced_cerebellum <- eBayes(fit.reduced_cerebellum, robust=TRUE)
summary(decideTests(fit.reduced_cerebellum))
diff_exp_cerebellum <-topTable(fit.reduced_cerebellum,coef=3,number=30000)

#let's subset - brain only
brain_only <- aggregated[grepl("brain",colnames(aggregated))]
brain_only_annot <- datSample[datSample$Tissue == "Brain",]
genotype_group_brain <- factor(brain_only_annot$Genotype,levels=c("WT","KO"))
sex_group_brain <- factor(brain_only_annot$Sex)
design_brain <- model.matrix(~ sex_group_brain + genotype_group_brain)
fit.reduced_brain <- lmFit(brain_only,design_brain)
fit.reduced_brain <- eBayes(fit.reduced_brain, robust=TRUE)
summary(decideTests(fit.reduced_brain))
diff_exp_brain <-topTable(fit.reduced_brain,coef=3,number=30000)

#let's subset - liver only
liver_only <- aggregated[grepl("liver",colnames(aggregated))]
liver_only_annot <- datSample[datSample$Tissue == "Liver",]
genotype_group_liver <- factor(liver_only_annot$Genotype,levels=c("WT","KO"))
sex_group_liver <- factor(liver_only_annot$Sex)
design_liver <- model.matrix(~sex_group_liver + genotype_group_liver )
fit.reduced_liver <- lmFit(liver_only,design_liver)
fit.reduced_liver <- eBayes(fit.reduced_liver, robust=TRUE)
summary(decideTests(fit.reduced_liver))
diff_exp_liver <-topTable(fit.reduced_liver,coef=3,number=30000)

#let's get lists of these
sig_dmrs_muscle <- data.frame(rownames(diff_exp_muscle[diff_exp_muscle$adj.P.Val < .01 &
                                                         abs(diff_exp_muscle$logFC) >.10,]))
sig_dmrs_heart <- data.frame(rownames(diff_exp_heart[diff_exp_heart$adj.P.Val < .01 &
                                                       abs(diff_exp_heart$logFC) >.10,]))
sig_dmrs_cerebellum <- data.frame(rownames(diff_exp_cerebellum[diff_exp_cerebellum$adj.P.Val < .01 &
                                                                 abs(diff_exp_cerebellum$logFC) >.10,]))
sig_dmrs_brain <- data.frame(rownames(diff_exp_brain[diff_exp_brain$adj.P.Val < .01 &
                                                       abs(diff_exp_brain$logFC) >.10,]))
sig_dmrs_liver <- data.frame(rownames(diff_exp_liver[diff_exp_liver$adj.P.Val < .01 &
                                                       abs(diff_exp_liver$logFC) >.10,]))

write.csv(sig_dmrs_muscle, "sig_dmrs_muscle.csv")
write.csv(sig_dmrs_heart, "sig_dmrs_heart.csv")
write.csv(sig_dmrs_cerebellum, "sig_dmrs_cerebellum.csv")
write.csv(sig_dmrs_brain, "sig_dmrs_brain.csv")
write.csv(sig_dmrs_liver, "sig_dmrs_liver.csv")

#Let's try to do a Venn Diagram.
a <- list('Muscle' = sig_dmrs_muscle[,1],
          'Heart' = sig_dmrs_heart[,1],
          'Cerebellum' = sig_dmrs_cerebellum[,1],
          'Brain' = sig_dmrs_brain[,1],
          'Liver' = sig_dmrs_liver[,1])
ggvenn(a,c("Muscle","Liver","Heart")) 

select(rownames(sig_dmrs_heart) %in% )

clock_results <- read_csv("data/ResultsClock/OutputMouseClockFeb8.2021.csv")

#Clock analysis looking at different tissues.
horvathclock_heart <- clock_results[clock_results$Tissue=="Heart",]
summarized_heart <- horvathclock_heart %>% group_by(Genotype) %>% summarise(HorvathClockAge = 
                                                                              mean(DNAmAgeFinalClockHeart), 
                                                        se= std.error(DNAmAgeFinalClockHeart))
summarized_heart$tissue <- "Heart"
horvathclock_muscle <- clock_results[clock_results$Tissue=="Muscle",]
summarized_muscle <- horvathclock_muscle %>% group_by(Genotype) %>% summarise(HorvathClockAge = 
                                                                                mean(DNAmAgeFinalClockMuscle), 
                                                                            se= std.error(DNAmAgeFinalClockMuscle))
summarized_muscle$tissue <- "Muscle"
horvathclock_liver <- clock_results[clock_results$Tissue=="Liver",]
summarized_liver <- horvathclock_liver %>% group_by(Genotype) %>% summarise(HorvathClockAge = 
                                                                              mean(DNAmAgeFinalClockLiver), 
                                                                             se= std.error(DNAmAgeFinalClockLiver))
summarized_liver$tissue <- "Liver"
horvathclock_brain <- clock_results[clock_results$Tissue=="Brain",]
summarized_brain <- horvathclock_brain %>% group_by(Genotype) %>% summarise(HorvathClockAge = 
                                                                              mean(DNAmAgeFinalClockBrain), 
                                                                            se= std.error(DNAmAgeFinalClockBrain))
summarized_brain$tissue <- "Brain"
horvathclock_cerebellum <- clock_results[clock_results$Tissue=="Cerebellum",]
summarized_cerebellum <- horvathclock_cerebellum %>% group_by(Genotype) %>% summarise(HorvathClockAge = 
                                                                                        mean(DNAmAgeFinalClockBrain), 
                                                                            se= std.error(DNAmAgeFinalClockBrain))
summarized_cerebellum$tissue <- "Cerebellum"
total_data <- rbind(summarized_heart,summarized_muscle,summarized_liver,summarized_brain,
                    summarized_cerebellum)
total_data$se <- total_data$se*100
total_data$HorvathClockAge <- total_data$HorvathClockAge*100
ggplot(total_data, aes(x=tissue,fill=Genotype,y=HorvathClockAge, 
                       ymin=HorvathClockAge-se, ymax = HorvathClockAge+se)) + 
  geom_bar(stat="identity",position=position_dodge(), color="black") + theme_minimal() +
  geom_errorbar(width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values=c('#FF9999','#FF0055')) + labs(x="Tissue", 
                                                          y = "Predicted Horvath Clock Age (% of lifespan)") +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"))
  
#Let's look at GO terms enriched for heart tissues only.
#GO analysis
sig_dmrs_heart <- data.frame(gsub("\\..*","",rownames(diff_exp_heart[diff_exp_heart$adj.P.Val <= .001 &
                                                                       abs(diff_exp_heart$logFC) >=.025,])))
not_sig_dmrs_heart <- data.frame(gsub("\\..*","",rownames(diff_exp_heart[diff_exp_heart$adj.P.Val > .001 &
                                                                           abs(diff_exp_heart$logFC) <.025,])))
colnames(sig_dmrs_heart) <- "gene_ids"
colnames(not_sig_dmrs_heart) <- "gene_ids"
sig_dmrs_heart$is_candidate <- 1
not_sig_dmrs_heart$is_candidate <- 0
combined_heart <- rbind(sig_dmrs_heart,not_sig_dmrs_heart)
res_hyper_bg_heart = go_enrich(combined_heart, n_randsets=100,organism="Mus.musculus")
top_gos_hyper_bg_heart = res_hyper_bg_heart[[1]][1:5, 'node_id']
plot_stats_heart = plot_anno_scores(res_hyper_bg_heart, top_gos_hyper_bg_heart)
plot_stats_heart
res_hyper_bg_heart[[1]][1:10, 'node_name']

#Volcano plot for heart only
cols <- densCols(diff_exp_heart$logFC, -log10(diff_exp_heart$adj.P.Val),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= diff_exp_heart$logFC, 
     y = -log10(diff_exp_heart$adj.P.Val), 
     col=cols, panel.first=grid(),
     main="Volcano plot of KO vs. WT mice DNA methylation levels at a CpG site (heart only)", 
     xlim=c(-.5,.5),
     ylim=c(0,10),
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     cex=0.6)
gn.selected <- abs(diff_exp_heart$logFC) >.16 & (diff_exp_heart$adj.P.Val) < .00001
text(diff_exp_heart$logFC[gn.selected],
     -log10(diff_exp_heart$adj.P.Val)[gn.selected],
     lab=rownames(diff_exp_heart)[gn.selected ], cex=0.6)

#Printing top differentially methylated CpGs for heart and overall.
heart_hits <- data.frame((diff_exp_heart[diff_exp_heart$adj.P.Val < .05,]))
hits <- data.frame(diff_exp[diff_exp$adj.P.Val < .05,])
write.csv(heart_hits,"diff_methexp_heart.csv")
write.csv(hits,"diff_methexp.csv")

heart <- orig[grepl("heart",colnames(orig))]
liver <- orig[grepl("liver",colnames(orig))]
heart_and_liver <- cbind(heart,liver)
heart_and_liver_annot <- datSample[datSample$Tissue == "Liver" | datSample$Tissue == "Heart",]
matched_positions <- match(rownames(heart_and_liver),annotations$CGid) #Finds valid matches in table.

#Looking at how CpGs at certain locations change.
matched_category <- annotations$main_Categories[matched_positions]
heart_and_liver <- drop_na(heart_and_liver)
group <- aggregate(x = heart_and_liver,                # Specify data column
                        by = list(matched_category),              # Specify group indicator
                        FUN = mean)
rotated_group <- data.frame(t(group))
colnames(rotated_group) <- rotated_group[1,]
rotated_group <- rotated_group[c(-1),]
melted_groups <- melt(as.matrix(rotated_group))
colnames(melted_groups) <- c("id","region","perc_meth")
melted_groups$KO <- grepl("KO",melted_groups$id)
melted_groups$Tissue <- substr(melted_groups$id,1,5)
melted_groups$perc_meth <- as.numeric(melted_groups$perc_meth)
summary <- melted_groups %>% group_by(region,KO,Tissue) %>% summarize(mean_perc = mean(perc_meth),
                                                                  sd_perc = sd(perc_meth))
colnames(summary) <- c("Region","KO","Tissue","Methylation_Perc","Methylation_SD")
summary$Region <- c(rep("Exon",4),rep("5' UTR",4),rep("IG Downstream",4),rep("IG Upstream",4),
                    rep("Intron",4),rep("Promoter",4),rep("3' UTR",4))
ggplot(summary,aes(x=Region,y=Methylation_Perc,fill=KO)) +
  geom_bar(stat="identity",position="dodge") + theme_classic() + 
  geom_errorbar(aes(ymin=Methylation_Perc-Methylation_SD,ymax=Methylation_Perc+Methylation_SD),width=.2,
                position=position_dodge(.9)) +
  facet_wrap(. ~ Tissue) +
  scale_fill_manual(values=c('lavender','red'))


