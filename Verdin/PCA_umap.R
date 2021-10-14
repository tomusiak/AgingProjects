library(readr)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(dplyr)
library(umap)
library(janitor)
library(corrplot)
library(reshape)
library(Hmisc)
library(tidyverse)
library(RColorBrewer)

setwd("/home/atom/Desktop/AgingProjects/Verdin/")
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
colnames(flipped_data)
ggplot(flipped_data, aes(age, fdg_high_temra_cd4s)) +
  geom_point() +
  geom_smooth(method=lm)

ggscatter(flipped_data, x = "age", y = "FDG_Pos_CM_CD8s", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Age", ylab = "% FDG High Central Memory CD8s ")

aged_data <- subset(flipped_data, select = -c(age) )
aged_data <- data.frame(scale(aged_data))
pca <-prcomp(aged_data)
fviz_eig(pca)
pca_plus_age <- cbind(pca,aged_data$age)
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             select.var = list(name =NULL, cos2 = NULL, contrib = 20)
  )
factors <- data.frame(pca$x)
factors$age <- flipped_data$age
ggscatter(factors, x = "age", y = "PC1",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Age", ylab = "PC1")

aged_umap <- umap(aged_data,labels = factors$age)
df <- data.frame(PC1 = aged_umap$layout[,1],
                 PC2 = aged_umap$layout[,2],
                 Age = factors$age)
ggscatter(df, x= "PC1", y="PC2", color="Age")
ggscatter(df, x = "Age", y = "PC2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Age", ylab = "PC2")

correlation <- cor(flipped_data)
correlation_test <- rcorr(as.matrix(flipped_data))
heatmap(correlation)
age_correlation <- data.frame(data.frame(correlation)[,1])
rownames(age_correlation) <- rownames(correlation)
age_correlation_test <- data.frame(correlation_test[3])
age_correlation_test <- data.frame(age_correlation_test[,1])
rownames(age_correlation_test) <- rownames(correlation)                                   
age_correlation <- tail(age_correlation,20)
age_correlation$Marker <- rownames(age_correlation)
age_correlation_test$Marker <- rownames(age_correlation_test)
age_correlation_test <- tail(age_correlation_test,20)
age_correlation$Significance <- age_correlation_test[,1]
colnames(age_correlation) <- c("Correlation", "Marker", "Significance")
age_correlation<- age_correlation %>% arrange(Correlation)
age_correlation$Significance <- round(age_correlation$Significance,3)
age_correlation$Significance <- paste("(", age_correlation$Significance,sep = "")
age_correlation$Significance <- paste(age_correlation$Significance,")",sep = "")
age_correlation$Correlation <- round(age_correlation$Correlation,2)
age_correlation$Combined <- paste(age_correlation$Correlation, age_correlation$Significance)

ggplot(age_correlation, aes(x = Correlation, y = reorder(Marker, -Correlation), fill=Correlation)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(label=age_correlation$Combined,hjust=-.02) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) + ylab("Marker") + xlab("Correlation") + guides(fill=FALSE) +
  xlim(c(0, .7)) + scale_fill_distiller(palette = "Reds", direction=11)

sup_learning <- read_csv("correlation_values.csv")
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE) / 4.47
      )
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
sup_learning_sum <- data_summary(sup_learning, varname="Data", 
                    groupnames=c("Model"))
ggplot(sup_learning_sum, aes(x=Model, y=Data, fill = "Red")) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=Data-se, ymax=Data+se), width=.2,
                position=position_dodge(.9)) + theme_classic() + scale_fill_brewer(palette = "Reds", direction=11) +
  ylab("Correlation with Age") + xlab("Model") + guides(fill=FALSE) + theme(axis.text=element_text(size=15))

