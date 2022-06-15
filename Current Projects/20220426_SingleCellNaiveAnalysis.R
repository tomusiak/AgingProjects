#Code to analyze some previously-acquired PBMC single-cell data.
#Goal is to focus on naive T cells and determine if 'naive T cell subsets' exist.
#If they exist, what markers do they uniquely express?

library(Seurat)
library(Matrix)
library(tidyverse)
library(GOfuncR)
library('biomaRt')
# 
# matrix <- readMM("Data/Witkowski2020/matrix.mtx")
# genes <- read_tsv("Data/Witkowski2020/genes.tsv", col_names = FALSE)$X1
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
#                                                           "hgnc_symbol"),values=genes,mart= mart)

pbmc <- readRDS("pbmc.rds")
True_Naive_T_Cells <- readRDS("truenaive.rds")

# match <- match(genes,G_list$ensembl_gene_id, nomatch=NULL)
# symbols <- G_list$hgnc_symbol[match]
# symbols <- make.unique(symbols)
# symbols[4] <- ".."
# symbols[6] <- "..."
# cell_ids <- read_tsv("Data/Witkowski2020/barcodes.tsv", col_names = FALSE)$X1
# rownames(matrix) <- symbols
# colnames(matrix) <- cell_ids
# pbmc <- CreateSeuratObject(counts = matrix)
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# top10 <- head(VariableFeatures(pbmc), 10)
# plot1 <- VariableFeaturePlot(pbmc)
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc)
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# DimPlot(pbmc, reduction = "pca")
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# ElbowPlot(pbmc)
# pbmc <- FindNeighbors(pbmc, dims = 1:10)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# pbmc <- RunUMAP(pbmc, dims = 1:10)
# saveRDS(pbmc, file = "pbmc.rds")

DimPlot(pbmc, reduction = "umap")
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# pbmc.markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)
# FeaturePlot(pbmc, features = c("CD3E", "CD8A", "CD4", "FOXP3", "CCR7", "GZMB", "PDCD1", "IL7R",
#                                "CD27"))
# Naive_T_Cells <- subset(pbmc, subset=seurat_clusters %in% c(0,3,5))
# DimPlot(Naive_T_Cells, reduction = "umap")
# FeaturePlot(Naive_T_Cells, features = c("CD8B", "IL2RG", "CD3E", "GZMA", "CD4", "IL7R", "FOXP3", "CCR7",
#                                "CD8A"))
# Naive_T_Cells <- NormalizeData(Naive_T_Cells, normalization.method = "LogNormalize", scale.factor = 10000)
# FindVariableFeatures(Naive_T_Cells, selection.method = "vst", nfeatures = 2000)
# Naive_T_Cells <- ScaleData(Naive_T_Cells)
# Naive_T_Cells <- RunPCA(Naive_T_Cells, features = VariableFeatures(object = Naive_T_Cells))
# DimPlot(Naive_T_Cells, reduction = "pca")
# DimHeatmap(Naive_T_Cells, dims = 1, cells = 500, balanced = TRUE)
# Naive_T_Cells <- FindNeighbors(Naive_T_Cells, dims = 1:10)
# Naive_T_Cells <- FindClusters(Naive_T_Cells, resolution = 0.3)
# Naive_T_Cells <- RunUMAP(Naive_T_Cells, dims = 1:10)
# DimPlot(Naive_T_Cells, reduction = "umap")
# FeaturePlot(Naive_T_Cells, features = c("CD8B", "IL32", "CD3E", "S100A4", "CD4", "JUNB", "KLF6", "CCR7",
#                                         "CD8A"))
# top10 <- head(VariableFeatures(Naive_T_Cells), 20)
# FindMarkers(Naive_T_Cells,ident.1=0, ident.2=2)

# True_Naive_T_Cells <- subset(Naive_T_Cells, subset=seurat_clusters %in% c(0,1))
# DimPlot(True_Naive_T_Cells, reduction = "umap")
# FeaturePlot(True_Naive_T_Cells, features = c("CD8B", "IL2RG", "CD3E", "GZMA", "CD4", "IL7R", "PTPRC", "CCR7",
#                                         "CD8A"))
# True_Naive_T_Cells <- NormalizeData(True_Naive_T_Cells, normalization.method = "LogNormalize", scale.factor = 10000)
# FindVariableFeatures(True_Naive_T_Cells, selection.method = "vst", nfeatures = 2000)
# True_Naive_T_Cells <- ScaleData(True_Naive_T_Cells)
# True_Naive_T_Cells <- RunPCA(True_Naive_T_Cells, features = VariableFeatures(object = True_Naive_T_Cells))
# DimPlot(True_Naive_T_Cells, reduction = "pca")
# DimHeatmap(True_Naive_T_Cells, dims = 1, cells = 500, balanced = TRUE)
# True_Naive_T_Cells <- FindNeighbors(True_Naive_T_Cells, dims = 1:10)
# True_Naive_T_Cells <- FindClusters(True_Naive_T_Cells, resolution = .6)
# True_Naive_T_Cells <- RunUMAP(True_Naive_T_Cells, dims = 1:10)
DimPlot(True_Naive_T_Cells, reduction = "umap",pt.size=.8)
FeaturePlot(True_Naive_T_Cells, features = c("CD3E", "CD4", "CCR7", "PDCD1", "SOX4", "ARPC1B","TMSB10","TOX","CD14"))
top50 <- head(VariableFeatures(True_Naive_T_Cells), 50)
#saveRDS(True_Naive_T_Cells, file = "truenaive.rds")
markers <- FindMarkers(True_Naive_T_Cells,ident.1=1, ident.2=2)
true_markers <- data.frame(rownames(markers)[markers$p_val_adj<.05])
gene_list <- data.frame(rownames(GetAssayData(object = pbmc, slot = "counts")))
true_markers$candidate <- 1
gene_list$candidate <- 0
colnames(true_markers) <- c("name","candidate")
colnames(gene_list) <- c("name","candidate")
full_list <- rbind(true_markers,gene_list)
go_results <- go_enrich(full_list)
go_results$results
write.csv(markers,"naive_scrnaseq_markers.csv")
                        
DoHeatmap(subset(True_Naive_T_Cells, downsample = 100), features = c("CD8B", "IL2RG", "CD3E", "GZMA", "CD4", "IL7R", "PTPRC", "CCR7", "CD8A"), size = 3)

                                                                     