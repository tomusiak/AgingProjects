library(readr)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(dplyr)
library(umap)
library(janitor)
setwd("/home/atom/Desktop/AgingProjects/")
data <- read_csv("20donors.csv", col_names=FALSE)
data <- data.frame(data)
new_data <- data.frame(data[,-c(1,3,4)])
new_data$X2[1] <- "age"
flipped_data <- data.frame(t(new_data))
colnames(flipped_data) <- unlist(flipped_data[1,])
flipped_data <- flipped_data[-1,]
flipped_data[ , ] <- apply(flipped_data[ , ], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
rownames(flipped_data) <- c(1:30)
flipped_data <- clean_names(flipped_data)
colnames(flipped_data)
ggplot(flipped_data, aes(age, live_t_cells_cd8_ro_neg_terminal_effectors_fdg_pos)) +
  geom_point() +
  geom_smooth(method=lm)

ggscatter(flipped_data, x = "age", y = "live_t_cells_cd8_ro_neg_naive_and_early_memory_early_memory_fdg_neg", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Age", ylab = "% SA-b-Gal Positive Effector Memory CD8s ")

aged_data <- subset(flipped_data, select = -c(age) )
aged_data <- data.frame(scale(aged_data))
pca <-prcomp(aged_data)
fviz_eig(pca)
pca_plus_age <- cbind(pca,aged_data$age)
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
factors <- data.frame(pca$x)
factors$age <- flipped_data$age
ggscatter(factors, x = "age", y = "PC1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Age", ylab = "PC1")

fdg_data <- select(aged_data,contains("FDG"))
fdg_pca <-prcomp(fdg_data)
fviz_eig(fdg_pca)
fviz_pca_ind(fdg_pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(fdg_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
factors <- data.frame(fdg_pca$x)
factors$age <- flipped_data$age
factors$nemfdgpos <- flipped_data$live_t_cells_cd8_ro_neg_naive_and_early_memory_early_memory_fdg_pos
ggscatter(factors, x = "PC1", y = "nemfdgpos", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PC1", ylab = "Age")
ggscatter(factors, x = "PC1", y = "PC2", color="nemfdgneg", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PC2", ylab = "Age")

aged_umap <- umap(aged_data,labels = factors$age)
df <- data.frame(PC1 = aged_umap$layout[,1],
                 PC2 = aged_umap$layout[,2],
                 Age = factors$age)
ggscatter(df, x = "PC2", y = "Age", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PC2", ylab = "Age")
