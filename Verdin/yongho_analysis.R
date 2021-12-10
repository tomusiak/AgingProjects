library(Rsamtools)
library(Rsubread)
library(DESeq2)
library(apeglm)
library(magrittr)
library(M3C)
library(readr)
library(dplyr)
library(umap)
library(stringr)
library(biomaRt)
library(tidyr)
library(ggplot2)
library(ashr)

setwd("/home/atom/Desktop/Data/Yong-Ho/BAMs")
BAMs <- list.files(pattern = "\\.BAM$")
#raw_counts <- featureCounts(BAMs,annot.inbuilt="hg38",isPairedEnd=TRUE,nthreads=5)
#raw_count_matrix <- raw_counts$counts
#write.csv(raw_count_matrix,"yongho_counts_raw.csv")
raw_count_matrix <- read.csv("~/Desktop/Data/Yong-Ho/BAMs/yongho_counts_raw.csv", row.names=1)
keep <- rowSums((raw_count_matrix)) >= 15
raw_count_matrix <- raw_count_matrix[keep,]
yongho_metadata <- read_csv("~/Desktop/Data/Yong-Ho/yongho_metadata.csv")
yongho_metadata <- yongho_metadata %>% arrange(id)
colnames(raw_count_matrix) <- yongho_metadata$id

gene_list <- rownames(raw_count_matrix)
#httr::set_config(httr::config(ssl_verifypeer = FALSE))
#id_to_name_mapping <- getBM(attributes = c("entrezgene_id","hgnc_symbol"), 
#                            filters = "entrezgene_id", values = gene_list,
#                            mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
#write.csv(id_to_name_mapping,"id_to_name_mapping.csv")
id_to_name_mapping <- read.csv("~/Desktop/Data/Yong-Ho/BAMs/id_to_name_mapping.csv", row.names=1)
matrix_ids <- rownames(raw_count_matrix)
matched_positions <- match(matrix_ids,id_to_name_mapping$entrezgene_id)
matched_symbols <- id_to_name_mapping$hgnc_symbol[matched_positions]
unmatched_positions <- match(matched_symbols,NA)
unmatched_positions_2 <- match(matched_symbols,"")
unmatched_positions <- unmatched_positions | unmatched_positions_2
unmatched_positions <- !is.na(unmatched_positions )
matched_symbols[unmatched_positions] <- matrix_ids[unmatched_positions]
unique_symbols <- make.names(matched_symbols, unique = TRUE)
rownames(raw_count_matrix) <- unique_symbols

dds <- DESeqDataSetFromMatrix(raw_count_matrix,
                                   colData = yongho_metadata,
                                   design = ~ young_or_aged + cd8_or_cd4 + cm_or_te + sabgal_high_or_low +
                                     sex + age + pooled + sabgal_high_or_low:cm_or_te)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cd8_or_cd4","cd8","cd4"))
resLFC <- lfcShrink(dds, res = res, contrast = "cd8_or_cd4",type="ashr")
plotMA(resLFC, ylim=c(-2,2))
resOrdered <- res[order(res$padj),]
top_hits_cd4 <- head(resOrdered,25)


res <- results(dds, contrast=c("sabgal_high_or_low","high","low"))
resLFC <- lfcShrink(dds, res = res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(resLFC, ylim=c(-2,2))
resOrdered <- resLFC[order(resLFC$padj),]
top_hits_sabgal <- head(resOrdered,1000)
write.csv(top_hits_sabgal,"yongho_top_hits_sabgal.csv")
top_hits_sabgal
head(top_hits_sabgal,10)


dds_norm <- vst(dds)
colnames(dds_norm)
normalized_counts <- assay(dds_norm) %>% t()
umap_results <- umap::umap(normalized_counts)
complete_umap <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("id") %>%
  dplyr::inner_join(yongho_metadata, by = "id")
ggplot(
  complete_umap,
  aes( x = X1, y = X2, color=cm_or_te)) +
  geom_point() + theme_classic() + labs(title="RNA-Seq UMAP Dimensional Reduction - Cell Type", x="Component 1", y="Component 2",
                                        color = "Central Memory or \n Terminal/Effector \n Memory") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1.5, 'cm'), legend.title = element_text(size=10))
ggplot(
  complete_umap,
  aes( x = X1, y = X2, color=donor_id_long)) +
  geom_point() + theme_classic() + labs(title="RNA-Seq UMAP Dimensional Reduction - Donors", x="Component 1", y="Component 2",
                                        color = "Donor ID                   ") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1.5, 'cm'), legend.title = element_text(size=10))

cm_samples <- colnames(dds_norm)[grepl("CM",colnames(dds_norm),fixed=TRUE)]
te_samples <- colnames(dds_norm)[grepl("TE",colnames(dds_norm),fixed=TRUE)]
cm_counts <- assay(dds_norm)[,cm_samples ] %>% t()
te_counts <-  assay(dds_norm)[,te_samples] %>% t()

umap_results_cm <- umap::umap(cm_counts)
cm_umap <- data.frame(umap_results_cm$layout) %>%
  tibble::rownames_to_column("id") %>%
  dplyr::inner_join(yongho_metadata, by = "id")
ggplot(
  cm_umap,
  aes( x = X1, y = X2, color=young_or_aged)) +
  geom_point() + theme_classic() + labs(title="RNA-Seq UMAP Dimensional Reduction CM Only - SAbGal", x="Component 1", y="Component 2",
                                        color = "SAbGal High Or Low") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1.5, 'cm'), legend.title = element_text(size=10))

umap_results_te <- umap::umap(te_counts)
te_umap <- data.frame(umap_results_te$layout) %>%
  tibble::rownames_to_column("id") %>%
  dplyr::inner_join(yongho_metadata, by = "id")
ggplot(
  te_umap,
  aes( x = X1, y = X2, color=young_or_aged)) +
  geom_point() + theme_classic() + labs(title="RNA-Seq UMAP Dimensional Reduction EM Only - SAbGal", x="Component 1", y="Component 2",
                                        color = "SAbGal High Or Low") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1.5, 'cm'), legend.title = element_text(size=10))

CD4_samples <- colnames(dds_norm)[grepl("CD4",colnames(dds_norm),fixed=TRUE)]
CD8_samples <- colnames(dds_norm)[grepl("CD8",colnames(dds_norm),fixed=TRUE)]
CD4_counts <- assay(dds_norm)[,CD4_samples ] %>% t()
CD8_counts <-  assay(dds_norm)[,CD8_samples] %>% t()

umap_results_CD4 <- umap::umap(CD4_counts)
CD4_umap <- data.frame(umap_results_CD4$layout) %>%
  tibble::rownames_to_column("id") %>%
  dplyr::inner_join(yongho_metadata, by = "id")
ggplot(
  CD4_umap,
  aes( x = X1, y = X2, color=sabgal_high_or_low)) +
  geom_point() + theme_classic() + labs(title="RNA-Seq UMAP Dimensional Reduction CD4 Only - SAbGal", x="Component 1", y="Component 2",
                                        color = "SAbGal High Or Low") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1.5, 'cm'), legend.title = element_text(size=10))

umap_results_CD8 <- umap::umap(CD8_counts)
CD8_umap <- data.frame(umap_results_CD8$layout) %>%
  tibble::rownames_to_column("id") %>%
  dplyr::inner_join(yongho_metadata, by = "id")
ggplot(
  CD8_umap,
  aes( x = X1, y = X2, color=sabgal_high_or_low)) +
  geom_point() + theme_classic() + labs(title="RNA-Seq UMAP Dimensional Reduction CD8 Only - SAbGal", x="Component 1", y="Component 2",
                                        color = "SAbGal High Or Low") +
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1.5, 'cm'), legend.title = element_text(size=10))

res <- results(dds)
resOrdered <- res[order(res$padj),]
top_hits_combo <- head(resOrdered,10)
write.csv(top_hits_sabgal,"top_hits_sabgal.csv")
pd1_with_age <- plotCounts(dds, gene="FOXP3", intgroup=c("cd8_or_cd4","sabgal_high_or_low"),returnData=TRUE)
ggplot(pd1_with_age, aes(x=sabgal_high_or_low, y=count,fill=cd8_or_cd4)) + 
  geom_violin() + theme_classic() + labs(title="CD57 Expression On T Cell Subsets in SA-bGal High Or Low Cells",x = "SAbGal High or Low",
  y="Counts",fill="Central Memory (CM) \n or Terminal + \n Effector Memory (TE)")

res <- results(dds, contrast=c("young_or_aged","young","aged"))
resLFC <- lfcShrink(dds, res = res, contrast = "young_or_aged",type="ashr")
plotMA(resLFC, ylim=c(-2,2))
resOrdered <- resLFC[order(res$padj),]
top_hits_aged <- head(resOrdered,1000)
write.csv(top_hits_aged,"top_hits_aged.csv")

raw_count_matrix["CD8A",]

CD4_samples <- colnames(raw_count_matrix)[grepl("CD4",colnames(dds_norm),fixed=TRUE)]
CD4_counts <- raw_count_matrix[,CD4_samples ]
rownames(yongho_metadata) <- yongho_metadata$id
CD4_metadata <-  yongho_metadata[CD4_samples,]
CD4_dds <- DESeqDataSetFromMatrix(CD4_counts,
                              colData = CD4_metadata,
                              design = ~ young_or_aged + cm_or_te + sabgal_high_or_low +
                                sex + age + pooled + sabgal_high_or_low:cm_or_te)
CD4_dds <- DESeq(CD4_dds)
CD4_res <- results(CD4_dds, contrast=c("sabgal_high_or_low","high","low"))
CD4_resLFC <- lfcShrink(CD4_dds, res = CD4_res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(CD4_resLFC, ylim=c(-2,2))
CD4_resOrdered <- CD4_resLFC[order(CD4_resLFC$padj),]
top_hits_cd4_sabgal <- head(CD4_resOrdered,250)

CD8_samples <- colnames(raw_count_matrix)[grepl("CD8",colnames(dds_norm),fixed=TRUE)]
CD8_counts <- raw_count_matrix[,CD8_samples ]
CD8_metadata <-  yongho_metadata[CD8_samples,]
CD8_dds <- DESeqDataSetFromMatrix(CD8_counts,
                                  colData = CD8_metadata,
                                  design = ~ young_or_aged + cm_or_te + sabgal_high_or_low +
                                    sex + age + pooled + sabgal_high_or_low:cm_or_te)
CD8_dds <- DESeq(CD8_dds)
CD8_res <- results(CD8_dds, contrast=c("sabgal_high_or_low","high","low"))
CD8_resLFC <- lfcShrink(CD8_dds, res = CD8_res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(CD8_resLFC, ylim=c(-2,2))
CD8_resOrdered <- CD8_resLFC[order(CD8_resLFC$padj),]
top_hits_cd8_sabgal <- head(CD8_resOrdered,250)

CD4TE_samples <- colnames(raw_count_matrix)[grepl("CD4TE",colnames(dds_norm),fixed=TRUE)]
CD4TE_counts <- raw_count_matrix[,CD4TE_samples ]
CD4TE_metadata <-  yongho_metadata[CD4TE_samples,]
CD4TE_dds <- DESeqDataSetFromMatrix(CD4TE_counts,
                                  colData = CD4TE_metadata,
                                  design = ~ young_or_aged + sabgal_high_or_low +
                                    sex + age + pooled)
CD4TE_dds <- DESeq(CD4TE_dds)
CD4TE_res <- results(CD4TE_dds, contrast=c("sabgal_high_or_low","high","low"))
CD4TE_resLFC <- lfcShrink(CD4TE_dds, res = CD4TE_res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(CD4TE_resLFC, ylim=c(-2,2))
CD4TE_resOrdered <- CD4TE_resLFC[order(CD4TE_resLFC$padj),]
top_hits_CD4TE_sabgal <- head(CD4TE_resOrdered,250)
write.csv(top_hits_CD4TE_sabgal,"yongho_top_hits_CD4TE_sabgal.csv")

CD4CM_samples <- colnames(raw_count_matrix)[grepl("CD4CM",colnames(dds_norm),fixed=TRUE)]
CD4CM_counts <- raw_count_matrix[,CD4CM_samples ]
CD4CM_metadata <-  yongho_metadata[CD4CM_samples,]
CD4CM_dds <- DESeqDataSetFromMatrix(CD4CM_counts,
                                  colData = CD4CM_metadata,
                                  design = ~ young_or_aged + sabgal_high_or_low +
                                    sex + age + pooled)
CD4CM_dds <- DESeq(CD4CM_dds)
CD4CM_res <- results(CD4CM_dds, contrast=c("sabgal_high_or_low","high","low"))
CD4CM_resLFC <- lfcShrink(CD4CM_dds, res = CD4CM_res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(CD4CM_resLFC, ylim=c(-2,2))
CD4CM_resOrdered <- CD4CM_resLFC[order(CD4CM_resLFC$padj),]
top_hits_CD4CM_sabgal <- head(CD4CM_resOrdered,250)
write.csv(top_hits_CD4CM_sabgal,"yongho_top_hits_CD4CM_sabgal.csv")

CD8TE_samples <- colnames(raw_count_matrix)[grepl("CD8TE",colnames(dds_norm),fixed=TRUE)]
CD8TE_counts <- raw_count_matrix[,CD8TE_samples ]
CD8TE_metadata <-  yongho_metadata[CD8TE_samples,]
CD8TE_dds <- DESeqDataSetFromMatrix(CD8TE_counts,
                                  colData = CD8TE_metadata,
                                  design = ~ young_or_aged + sabgal_high_or_low +
                                    sex + age + pooled)
CD8TE_dds <- DESeq(CD8TE_dds)
CD8TE_res <- results(CD8TE_dds, contrast=c("sabgal_high_or_low","high","low"))
CD8TE_resLFC <- lfcShrink(CD8TE_dds, res = CD8TE_res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(CD8TE_resLFC, ylim=c(-2,2))
CD8TE_resOrdered <- CD8TE_resLFC[order(CD8TE_resLFC$padj),]
top_hits_CD8TE_sabgal <- head(CD8TE_resOrdered,250)
write.csv(top_hits_CD8TE_sabgal,"yongho_top_hits_CD8TE_sabgal.csv")

CD8CM_samples <- colnames(raw_count_matrix)[grepl("CD8CM",colnames(dds_norm),fixed=TRUE)]
CD8CM_counts <- raw_count_matrix[,CD8CM_samples ]
CD8CM_metadata <-  yongho_metadata[CD8CM_samples,]
CD8CM_dds <- DESeqDataSetFromMatrix(CD8CM_counts,
                                  colData = CD8CM_metadata,
                                  design = ~ young_or_aged + sabgal_high_or_low +
                                    sex + age + pooled)
CD8CM_dds <- DESeq(CD8CM_dds)
CD8CM_res <- results(CD8CM_dds, contrast=c("sabgal_high_or_low","high","low"))
CD8CM_resLFC <- lfcShrink(CD8CM_dds, res = CD8CM_res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(CD8CM_resLFC, ylim=c(-2,2))
CD8CM_resOrdered <- CD8CM_resLFC[order(CD8CM_resLFC$padj),]
top_hits_CD8CM_sabgal <- head(CD8CM_resOrdered,250)
write.csv(top_hits_CD8CM_sabgal,"yongho_top_hits_CD8CM_sabgal.csv")
