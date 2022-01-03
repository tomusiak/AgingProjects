#Sets location of Olga's data.
setwd("/home/atom/Desktop/Data/Olga") #Sets directory.

#Libraries to import.
library(WGCNA)
library(missMethyl)
library(limma)
library(dplyr)
require("biomaRt")
library("RColorBrewer")

#Helper function to go from human to mouse genes
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genesV2)
}

#Loading in data.
load("./data/NormalizedData/all_probes_sesame_normalized.Rdata")
datSample=read.csv("./data/SampleSheetAgeN80version4.csv")
annotations <- readRDS("./tools/geneAnnotations/AnnotationAminHaghani/Latest versions/Human.Homo_sapiens.hg38.Amin.V4.RDS", refhook = NULL)
dat0=data.frame(normalized_betas_sesame)

#Some QC checks and re-formatting performed by Horvath's group and copied here.
#QC primarily performed by cluster analysis.
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

#Setting gene names in place of CpG ID.
matched_positions <- match(dat0$CGid,annotations$CGid) #Finds valid matches in table.
matched_symbols <- annotations$SYMBOL[matched_positions]
dat0$symbol <- matched_symbols
unique_symbols <- make.names(matched_symbols, unique = TRUE)
rownames(dat0) <- unique_symbols
dat0 <- dat0[,2:81]

#limma differential analysis
genotype_group <- factor(datSample$Genotype,levels=c("WT","KO"))
tissue_group <- factor(datSample$Tissue)
design <- model.matrix(~genotype_group + tissue_group)
fit.reduced <- lmFit(dat0,design)
fit.reduced <- eBayes(fit.reduced, robust=TRUE)
summary(decideTests(fit.reduced))
top<-topTable(fit.reduced,coef=2,number=20000)
top_genes <- rownames(top)
top_genes <- gsub("\\..*","", top_genes )
matched_positions <- match(top_genes,genes$HGNC.symbol) #Finds valid matches in table.
matched_symbols <- genes$MGI.symbol[matched_positions]
rownames(top) <- make.names(matched_symbols,unique=TRUE)
cols <- densCols(top$logFC, -log10(top$adj.P.Val),
                 nbin=25, bandwidth=1,
                 colramp = colorRampPalette(brewer.pal(5, "Reds")))
plot(x= top$logFC, 
     y = -log10(top$adj.P.Val), 
     col=cols, panel.first=grid(),
     main="Volcano plot of KO vs. WT mice", 
     xlim=c(-.20,.20),
     ylim=c(0,20),
     xlab="Effect size: log(fold-change)",
     ylab="-log10(adjusted p-value)",
     cex=0.6)
gn.selected <- abs(top$logFC) >.025 & (top$adj.P.Val) < .0000005
text(top$logFC[gn.selected],
     -log10(top$adj.P.Val)[gn.selected],
     lab=rownames(top)[gn.selected ], cex=0.5)


mouse_genes <- matched_symbols[!is.na(matched_symbols)]
top_genes_unique <- (data.frame(top_genes) %>% distinct())[,1]
top["Kdm5a",]
rownames(dat0)
genes <- convertHumanGeneList(top_genes_unique)
matched_positions <- match(top_genes,genes$HGNC.symbol) #Finds valid matches in table.
matched_symbols <- genes$MGI.symbol[matched_positions]
mouse_genes <- matched_symbols[!is.na(matched_symbols)]
write.csv(mouse_genes,"mouse_genes.csv")

