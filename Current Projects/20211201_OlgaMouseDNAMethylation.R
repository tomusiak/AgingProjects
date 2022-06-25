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
library(stringr)
library(umap)
library(ggrepel)
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
# rownames(dat0) <- make.names(dat0$Symbol,unique=TRUE)
# dat0 <- dat0[,1:81]
# m_orig <- log2(orig/(1-orig))
# colnames(dat0)[-1]=paste0(datSample$OriginalOrderInBatch,datSample$Tissue)
# colnames(orig)=paste0(datSample$OriginalOrderInBatch,datSample$Tissue)

#processing
m_values <- log2(dat0/(1-dat0))

#Some QC checks and re-formatting performed by Horvath's group and copied here.
#QC primarily performed by cluster analysis
corSample=cor(dat0, use="p")
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

m_values$Symbol <- matched_symbols
m_values$category <- matched_category
m_values <- subset(m_values, m_values$category == "Exon" | m_values$category == "Promoter")
m_values <- m_values[!is.na(m_values$Symbol),]
dat0_symbols <- dat0
dat0_symbols$symbol <- matched_symbols
dat0_symbols$category <- matched_category
dat0_symbols <- subset(dat0_symbols, dat0_symbols$category == "Exon" | dat0_symbols$category == "Promoter")
dat0_symbols <- dat0_symbols[!is.na(dat0_symbols$symbol),]

#Let's do a complicated thing where we decipher most "significant" changes in
#CpGs per gene by doing column-wise subtraction. Heart only here.
genes_cpg_subtractions <- data.frame(matrix(ncol = 3, nrow = 1))
colnames(genes_cpg_subtractions) <- c("Gene","Mean","SE")
heart_m_values <- dat0_symbols[,grepl("heart",colnames(dat0_symbols)) |
                             grepl("symbol",colnames(dat0_symbols)) |
                               grepl("liver",colnames(dat0_symbols))]
split_m_values <- heart_m_values[,1:(ncol(heart_m_values)-1)] %>% 
  split(heart_m_values,f=heart_m_values$symbol)
transposed_m_values <- lapply(split_m_values,t)
dataframed_m_values <- lapply(transposed_m_values,as.data.frame)
for (i in 1:length(dataframed_m_values)) {
  dataframed_m_values[[i]]$Genotype <- str_sub(rownames((dataframed_m_values[[i]])),-2,) 
}
for (x in 1:length(dataframed_m_values)) {
  name <- names(dataframed_m_values[x])[[1]]
  practice_split<-split(dataframed_m_values[[x]],dataframed_m_values[[x]]$Genotype)
  for (i in 1:length(practice_split)) {
    practice_split[[i]] <- (practice_split[[i]])[,1:ncol((practice_split[[i]]))-1]
  }
  average_differences <- data.frame(practice_split[[2]] - practice_split[[1]])
  average_average_differences <- rowMeans(average_differences)
  average_average_sd <- sd(average_average_differences)
  average_average_average_differences <- mean(average_average_differences)
  average_average_n <- ncol(average_differences)
  if(average_average_n > 1) {
    genes_cpg_subtractions <- genes_cpg_subtractions %>% add_row(Gene = name, Mean = average_average_average_differences, 
                                     SE =average_average_sd/sqrt(average_average_n) )
  }
}
genes_cpg_subtractions <- genes_cpg_subtractions[-1,]
rownames(genes_cpg_subtractions) <- genes_cpg_subtractions[,1]
positives <- genes_cpg_subtractions[genes_cpg_subtractions$Mean > 0,]
ordered_list_positives <- positives[order(-(positives$Mean-abs(positives$SE))),]
negatives <- genes_cpg_subtractions[genes_cpg_subtractions$Mean < 0,]
ordered_list_negatives <- negatives[order((negatives$Mean+abs(negatives$SE))),]
top_changes <- rbind(ordered_list_positives[1:15,],ordered_list_negatives[1:15,])
ggplot(ordered_list_positives[1:15,],aes(x=reorder(Gene,-SE),y=Mean,ymin=Mean-SE,ymax=Mean+SE)) + 
  geom_bar(stat="identity",fill="maroon",color="red") + theme_classic() +
  geom_errorbar(width=.2,position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(ordered_list_negatives[1:15,],aes(x=reorder(Gene,-SE),y=Mean,ymin=Mean-SE,ymax=Mean+SE)) + 
  geom_bar(stat="identity",fill="maroon",color="red") + theme_classic() +
  geom_errorbar(width=.2,position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(top_changes,aes(x=reorder(Gene,-Mean),y=Mean,ymin=Mean-SE,ymax=Mean+SE)) + 
  geom_bar(stat="identity",fill="pink",color="red") + theme_classic() +
  geom_errorbar(width=.2,position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Gene") + ylab("Mean Change of Methylation (WT-KO)") +
  ggtitle("Mean Change in Methylation Site State, Per Gene (Heart and Liver Only)")
 
aggregated <- aggregate(x = m_values[,1:(ncol(m_values)-2)],                # Specify data column
          by = list(m_values$Symbol),              # Specify group indicator
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
qplot(c(1:20), var_explained[1:10]) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1) + theme_classic()

#PCA for heart and liver only.
rotated <- t(aggregated[grepl("heart",colnames(aggregated)) | grepl("liver",colnames(aggregated))])
pca<-prcomp(rotated,scale=TRUE)
var_explained = pca$sdev^2 / sum(pca$sdev^2)
pca <- data.frame(pca$x)
pca$genotype <- datSample[datSample$Tissue=="Heart" | datSample$Tissue=="Liver",]$Genotype
pca$tissue <- datSample[datSample$Tissue=="Heart"| datSample$Tissue=="Liver",]$Tissue
pca$sex <- datSample[datSample$Tissue=="Heart"| datSample$Tissue=="Liver",]$Sex
pca$age <- datSample[datSample$Tissue=="Heart"| datSample$Tissue=="Liver",]$Age
ggplot(data=data.frame(pca),aes(x=PC1,y=PC3,color=genotype,shape=tissue)) +
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

#Volcano plot colored by number of CpGs per gene.
symbols <- m_values$Symbol
symbol_occurrences <- data.frame(table(symbols))
rownames(symbol_occurrences) <- symbol_occurrences[,1]
symbol_occurrences <- subset(symbol_occurrences, select = c(Freq))
diff_exp_freq <- merge(diff_exp,symbol_occurrences,by="row.names")
rownames(diff_exp_freq) <- diff_exp_freq$Row.names

diff_exp_freq$label = abs(diff_exp_freq$logFC) >.06 & (diff_exp_freq$adj.P.Val) < .0000005

ggplot(diff_exp_freq,aes(x=logFC,y=-log10(adj.P.Val),size=sqrt(Freq),col=1/(sqrt(Freq)))) +
  geom_point() + theme_classic() + geom_text_repel(data = subset(diff_exp_freq, label==TRUE), aes(label = Row.names), 
                                                   box.padding = unit(0.60, "lines")) +
  scale_color_distiller(palette = "RdPu")

#Let's do the same but for heart and liver only.
heart_and_liver_aggregated <- aggregated[grepl("heart",colnames(aggregated)) | 
                                           grepl("liver",colnames(aggregated))]
rownames(heart_and_liver_aggregated) <- aggregated[,1]
heart_and_liver_annot <- datSample[datSample$Tissue == "Heart" | datSample$Tissue == "Liver",]
genotype_group_heartliver <- factor(heart_and_liver_annot$Genotype,levels=c("WT","KO"))
sex_group_heartliver <- factor(heart_and_liver_annot$Sex)
tissue_group_heartliver <- factor(heart_and_liver_annot$Tissue)
design_heartliver <- model.matrix(~ sex_group_heartliver + tissue_group_heartliver +
                                    genotype_group_heartliver)
fit.reduced_heartliver <- lmFit(heart_and_liver_aggregated,design_heartliver)
fit.reduced_heartliver <- eBayes(fit.reduced_heartliver, robust=TRUE)
summary(decideTests(fit.reduced_heartliver))
diff_exp_heartliver <-topTable(fit.reduced_heartliver,coef=3,number=30000)
diff_exp_heartliver_freq <- merge(diff_exp_heartliver,symbol_occurrences,by="row.names")
rownames(diff_exp_heartliver_freq) <- diff_exp_heartliver_freq$Row.names

diff_exp_heartliver_freq$label = abs(diff_exp_heartliver_freq$logFC) >2.5 & 
  (diff_exp_heartliver_freq$adj.P.Val) < .00000000000000000000001

ggplot(diff_exp_heartliver_freq,aes(x=logFC,y=-log10(adj.P.Val),size=sqrt(Freq),col=1/(sqrt(Freq)))) +
  geom_point() + theme_classic() + geom_text_repel(data = subset(diff_exp_heartliver_freq, label==TRUE), aes(label = Row.names), 
                                                   box.padding = unit(0.60, "lines")) +
  xlim(-5.5,5.5) + ylim(0,36) +
  scale_color_distiller(palette = "RdPu")

#Let's look at Slc10a1
Slc10a1_heartliver <- dat0_symbols[grepl("heart",colnames(dat0_symbols)) | 
                                           grepl("liver",colnames(dat0_symbols)) |
                                    grepl("symbol",colnames(dat0_symbols))]
Slc10a1_heartliver <- data.frame(t(Slc10a1_heartliver[Slc10a1_heartliver$symbol=="Slc10a1",
                                                  1:(ncol(Slc10a1_heartliver)-2)]))
Slc10a1_heartliver_melted <- melt(as.matrix(Slc10a1_heartliver))
Slc10a1_heartliver_melted$KO <- grepl("KO",Slc10a1_heartliver_melted$Var1)
Slc10a1_heartliver_summary <- Slc10a1_heartliver_melted %>% group_by(KO,Var2) %>% 
                                            summarize(mean_perc = mean(value),
                                                                    sd_perc = sd(value))
split_Slc10a1 <- split(Slc10a1_heartliver_summary,f=Slc10a1_heartliver_summary$KO)
cpg_change <- as.data.frame(split_Slc10a1[2])$TRUE.mean_perc -
  as.data.frame(split_Slc10a1[1])$FALSE.mean_perc 
ggplot(Slc10a1_heartliver_summary,aes(x=KO,y=mean_perc,group=Var2)) +
  geom_line() + geom_point() +
  theme_classic()

#Let's look at Dnajc22
Dnajc22_heartliver <- dat0_symbols[grepl("heart",colnames(dat0_symbols)) | 
                                    grepl("liver",colnames(dat0_symbols)) |
                                    grepl("symbol",colnames(dat0_symbols))]
Dnajc22_heartliver <- data.frame(t(Dnajc22_heartliver[Dnajc22_heartliver$symbol=="Dnajc22",
                                                    1:(ncol(Dnajc22_heartliver)-2)]))

Dnajc22_heartliver_melted <- melt(as.matrix(Dnajc22_heartliver))
Dnajc22_heartliver_melted$KO <- grepl("KO",Dnajc22_heartliver_melted$Var1)
Dnajc22_heartliver_summary <- Dnajc22_heartliver_melted %>% group_by(KO,Var2) %>% 
  summarize(mean_perc = mean(value),
            sd_perc = sd(value))
split_Dnajc22 <- split(Dnajc22_heartliver_summary,f=Dnajc22_heartliver_summary$KO)
cpg_change <- as.data.frame(split_Dnajc22[2])$TRUE.mean_perc -
  as.data.frame(split_Dnajc22[1])$FALSE.mean_perc 
ggplot(Dnajc22_heartliver_summary,aes(x=KO,y=mean_perc,group=Var2)) +
  geom_line() + geom_point() +
  theme_classic()


#Let's look at Id2
Id2_heartliver <- dat0_symbols[grepl("heart",colnames(dat0_symbols)) | 
                                     grepl("liver",colnames(dat0_symbols)) |
                                     grepl("symbol",colnames(dat0_symbols))]
Id2_heartliver <- data.frame(t(Id2_heartliver[Id2_heartliver$symbol=="Id2",
                                                      1:(ncol(Id2_heartliver)-2)]))

Id2_heartliver_melted <- melt(as.matrix(Id2_heartliver))
Id2_heartliver_melted$KO <- grepl("KO",Id2_heartliver_melted$Var1)
Id2_heartliver_summary <- Id2_heartliver_melted %>% group_by(KO,Var2) %>% 
  summarize(mean_perc = mean(value),
            sd_perc = sd(value))
split_Id2 <- split(Id2_heartliver_summary,f=Id2_heartliver_summary$KO)
cpg_change <- as.data.frame(split_Id2[2])$TRUE.mean_perc -
  as.data.frame(split_Id2[1])$FALSE.mean_perc 
ggplot(Id2_heartliver_summary,aes(x=KO,y=mean_perc,group=Var2)) +
  geom_line() + geom_point() +
  theme_classic()


