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
yongho_metadata <- read_csv("~/Desktop/Data/Yong-Ho/yongho_metadata.csv")
yongho_metadata <- yongho_metadata %>% arrange(id)
colnames(raw_count_matrix) <- yongho_metadata$id
filtered_count_matrix <- raw_count_matrix %>%
  dplyr::filter(rowSums(.) >= 15)

gene_list <- rownames(filtered_count_matrix)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
id_to_name_mapping <- getBM(attributes = c("entrezgene_id","hgnc_symbol"), 
                            filters = "entrezgene_id", values = gene_list,
                            mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
matrix_ids <- rownames(filtered_count_matrix)
matched_positions <- match(matrix_ids,id_to_name_mapping$entrezgene_id)
matched_symbols <- id_to_name_mapping$hgnc_symbol[matched_positions]
unmatched_positions <- match(matched_positions,NA)
matched_symbols[unmatched_positions] < matrix_ids[unmatched_positions]
unmatched_positions <- is.na(matched_symbols)
unmatched_symbols <- matrix_ids[unmatched_positions]
matched_symbols[unmatched_positions] <- unmatched_symbols
unique_symbols <- make.names(matched_symbols, unique = TRUE)
rownames(filtered_count_matrix) <- unique_symbols

dds <- DESeqDataSetFromMatrix(filtered_count_matrix,
                                   colData = yongho_metadata,
                                   design = ~ young_or_aged + cd8_or_cd4 + cm_or_te + sabgal_high_or_low +
                                     sex + age + pooled + young_or_aged:sabgal_high_or_low)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cd8_or_cd4","cd8","cd4"),addMLE=TRUE)
resOrdered <- res[order(res$padj),]
top_hits_cd4 <- head(resOrdered,25)

res <- results(dds, contrast=c("sabgal_high_or_low","high","low"))
resOrdered <- res[order(res$padj),]
top_hits_sabgal <- head(resOrdered,10)

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

length(te_counts)

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

res <- results(dds, contrast=c("sabgal_high_or_low","high","low"))
resLFC <- lfcShrink(dds, res = res, contrast = "sabgal_high_or_low",type="ashr")
plotMA(resLFC, ylim=c(-2,2))
resOrdered <- resLFC[order(res$padj),]
top_hits_sabgal <- head(resOrdered,100)
top_hits_sabgal

res <- results(dds)
resOrdered <- res[order(res$padj),]
top_hits_combo <- head(resOrdered,10)
write.csv(top_hits_sabgal,"top_hits_sabgal.csv")
pd1_with_age <- plotCounts(dds, gene="IL6", intgroup="sabgal_high_or_low",returnData=TRUE)
ggplot(pd1_with_age, aes(x=sabgal_high_or_low, y=count)) + 
  geom_bar(stat="identity", fill="gray") + theme_classic()

res <- results(dds, contrast=c("young_or_aged","young","aged"))
resLFC <- lfcShrink(dds, res = res, contrast = "young_or_aged",type="ashr")
plotMA(resLFC, ylim=c(-2,2))
resOrdered <- resLFC[order(res$padj),]
top_hits_aged <- head(resOrdered,100)
write.csv(top_hits_aged,"top_hits_aged.csv")
