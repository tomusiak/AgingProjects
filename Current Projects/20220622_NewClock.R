#This program creates a new epigenetic clock by training on CpGs that were identified to be not associated with
# T cell differentiation in a previous analysis. The clock is based on elastic net.

#Some useful code.
source("AgingProjects/Useful Scripts/generally_useful.R")
setwd("Data/") #Sets directory.

#Reading in necessary libraries.
library(glmnet)
library(readr)
library(tidyr)
library(caret)
library(MASS)

#Reading in database data for training and testing.
cpg_table <- read.csv("ClockConstruction/cpg_table.csv",row.names=1)
permitted_cpgs <- read.csv("ClockConstruction/nodiff_cpgs.csv",row.names=1)[,1]
sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
healthy_samples <- sample_table[sample_table$Condition=="Control",1]
healthy_sample_table <- sample_table[sample_table$Condition=="Control",]
cpg_table <- cpg_table[,colnames(cpg_table) %in% healthy_samples]
cpg_table_rotated <- data.frame(t(cpg_table))

#Annotating and splitting data.
cpg_table_rotated$Age <- healthy_sample_table$Age
cpg_table_rotated <- cpg_table_rotated[,colnames(cpg_table_rotated) %in%
                                  permitted_cpgs | colnames(cpg_table_rotated) ==
                                  "Age"]
splitting <- createDataPartition(cpg_table_rotated$Age,p=.67,list=FALSE)
training_set <- cpg_table_rotated[splitting,]
test_set <- cpg_table_rotated[-splitting,]

#Training the model.
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              search = "random",
                              verboseIter = TRUE)
model <-train(Age ~ .,
      data = training_set,
      method = "glmnet",
      preProcess = c("center", "scale"),
      tuneLength = 25,
      trControl = train_control)
predicted_age <- predict(model,test_set)
predicted_age <- data.frame(predicted_age)
predicted_age$ID <- rownames(predicted_age)
sample_table_2 <- merge(sample_table,predicted_age,by="ID")

#Useful script found online to create a ggplot with linear regression information.
ggplotRegression <- function(fit){
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    theme_classic() +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) 
}

#Plotting accuracy of the model.
ggplotRegression(lm(Age ~ predicted_age, sample_table_2))

#Identifying which CpGs are important
important_cpgs<-varImp(model)$importance
important_cpgs <- important_cpgs[order(important_cpgs$Overall,decreasing=TRUE),,drop=FALSE]

#Reading in validation table of data from Jonkman 2022. Will be useful for comparing clock performance.
validation_cpg_table <- read.csv("ClockConstruction/validation_cpg_table.csv",row.names=1)
validation_cpg_table <- validation_cpg_table[rownames(validation_cpg_table) %in% 
                                               rownames(cpg_table),]
validation_cpg_table_rotated <- data.frame(t(validation_cpg_table))
validation_cpg_table_rotated <- validation_cpg_table_rotated[,colnames(validation_cpg_table_rotated) %in%
                                                               permitted_cpgs]

#Fitting the Jonkman useful datasets to the model and exporting them to analyze whether they work for
# differentiation.
fitted_validation <- predict(model,validation_cpg_table_rotated)
write.csv(fitted_validation,"ClockConstruction/predictions.csv")

