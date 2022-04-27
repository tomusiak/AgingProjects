library(Seurat)
young_data <- data.frame(read_csv("~/Desktop/Data/GSE146395_young_plasma_count.csv"))
normal_data <- data.frame(read_csv("~/Desktop/Data/GSE146395_normal_aging_count.csv"))
aged_data <- data.frame(read_csv("~/Desktop/Data/GSE146395_aged_plasma_count.csv"))
rownames(young_data) <- young_data[,1]
rownames(normal_data) <- normal_data[,1]
rownames(aged_data) <- aged_data[,1]
#young_seurat <- NormalizeData(CreateSeuratObject(counts = young_data, project = "humanin", min.cells = 10, min.features = 400))
young_seurat <- AddMetaData(
  object = young_seurat,
  metadata = "young_plasma",
  col.name = 'age'
)
#normal_seurat <- NormalizeData(CreateSeuratObject(counts = normal_data, project = "humanin", min.cells = 10, min.features = 400))
normal_seurat <- AddMetaData(
  object = normal_seurat,
  metadata = "normal_plasma",
  col.name = 'age'
)
#aged_seurat <- NormalizeData(CreateSeuratObject(counts = aged_data, project = "humanin", min.cells = 10, min.features = 400))
aged_seurat <- AddMetaData(
  object = aged_seurat,
  metadata = "aged_plasma",
  col.name = 'age'
)
total_seurat <- merge(young_seurat, y = c(normal_seurat, aged_seurat), add.cell.ids = c("young_plasma", "normal_plasma", "aged_plasma"), project = "humanin", merge.data = FALSE)
total_seurat <- NormalizeData(total_seurat)
all.genes <- rownames(total_seurat)
total_seurat <- ScaleData(total_seurat, features = all.genes)
total_seurat <- FindVariableFeatures(total_seurat, selection.method = "vst", nfeatures = 2000)
total_seurat <- FindNeighbors(total_seurat, dims = 1:10)
total_seurat <- FindClusters(total_seurat, resolution = 0.5)
total_seurat <- RunPCA(total_seurat, features = VariableFeatures(object = total_seurat))
total_seurat <- RunUMAP(total_seurat, dims = 1:10)
DimPlot(total_seurat, reduction = "umap")
DimPlot(total_seurat, label=T, group.by="age", cells.highlight= list(young, old), cols.highlight = c("darkblue", "darkred"), cols= "grey")
Idents(total_seurat) <- "age"
mito_plot <- VlnPlot(object = total_seurat, features = c("mt-Rnr1","mt-Rnr2","Gm20594"))
inflam_plot <- VlnPlot(object = total_seurat, features = c("Junb","Icam1","Fos"))
markers <- FindMarkers(total_seurat, ident.1 = "young_plasma", ident.2 = "aged_plasma")


