#Sets location of Olga's data.
setwd("/home/atom/Desktop/Data/Olga") #Sets directory.

#Libraries to import.
library(WGCNA)
library(missMethyl)
library(limma)
library(dplyr)
require("biomaRt")
library("RColorBrewer")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(DMRcate)
library(readxl)

#Loading in data.
load("./data/NormalizedData/all_probes_sesame_normalized.Rdata")
datSample=read.csv("./data/SampleSheetAgeN80version4.csv")
annotations <-  readRDS("/home/atom/Desktop/Data/Olga/tools/geneAnnotations/AnnotationAminHaghani/Latest versions/Mouse.Mus_musculus.GRCm38.100.Amin.V6.RDS")
dat0=data.frame(normalized_betas_sesame)

#Let's swap CG ID's for gene names, and remove any NA values.
matched_positions <- match(dat0$CGid,annotations$CGid) #Finds valid matches in table.
matched_symbols <- annotations$SYMBOL[matched_positions]
dat0$Symbol <- matched_symbols
dat0 <- dat0[!is.na(dat0$Symbol),]
rownames(dat0) <- make.names(dat0$Symbol,unique=TRUE)
dat0 <- dat0[,1:81]

#Some QC checks and re-formatting performed by Horvath's group and copied here.
#QC primarily performed by cluster analysis
colnames(dat0)[-1]=paste0(datSample$OriginalOrderInBatch,datSample$Tissue)
corSample=cor(dat0[,-1], use="p")
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
mouse_list <- dat0[,2:81]
m_values <- log2(mouse_list/(1-mouse_list))

#limma differential methylation analysis
genotype_group <- factor(datSample$Genotype,levels=c("WT","KO"))
tissue_group <- factor(datSample$Tissue,levels=c("Heart","Brain","Liver","Muscle","Cerebellum"))
sex_group <- factor(datSample$Sex)
design <- model.matrix(~genotype_group + tissue_group + sex_group)
fit.reduced <- lmFit(m_values,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
diff_exp <-topTable(fit.reduced,coef=2,number=30000)

#Differential variability analysis using missMethyl.
fitvar <- varFit(m_values, design = design, coef = c(1,4))
summary(decideTests(fitvar))
topDV <- topVar(fitvar, coef=2,number=300)
diff_exp["Clk1.1",]

#let's pull only unique stuff
single_genes <- gsub("\\..*","", rownames(diff_exp) )
first_iteration <- diff_exp[!duplicated(single_genes,MARGIN=1),]

#limma differential analysis
cols <- densCols(first_iteration$logFC, -log10(first_iteration$adj.P.Val),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= first_iteration$logFC, 
     y = -log10(first_iteration$adj.P.Val), 
     col=cols, panel.first=grid(),
     main="Volcano plot of KO vs. WT mice", 
     xlim=c(-2,2),
     ylim=c(0,20),
     xlab="Effect size: log(fold-change)",
     ylab="-log10(adjusted p-value)",
     cex=0.6)
gn.selected <- abs(first_iteration$logFC) >.01 & (first_iteration$adj.P.Val) < .0000005
text(first_iteration$logFC[gn.selected],
     -log10(first_iteration$adj.P.Val)[gn.selected],
     lab=rownames(first_iteration)[gn.selected ], cex=0.6)

#other stuff
matched_positions <- match(single_genes,genes$HGNC.symbol) #Finds valid matches in table.
matched_symbols <- genes$MGI.symbol[matched_positions]
mouse_genes <- matched_symbols[!is.na(matched_symbols)]
matched_positions <- match(top_genes,genes$HGNC.symbol) #Finds valid matches in table.
matched_symbols <- genes$MGI.symbol[matched_positions]
rownames(top) <- make.names(make.unique(matched_symbols,sep="."))
single_genes <- gsub("\\..*","", rownames(top) )
first_iteration <- top[!duplicated(single_genes,MARGIN=1),]
