source("AgingProjects/Useful Scripts/generally_useful.R")
setwd("Data/") #Sets directory.

library(readr)
library(tidyr)
library(caret)
library(glmnet)
library(MASS)

sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
cpg_table <- read.csv("ClockConstruction/cpg_table.csv",row.names=1)
cpg_table_rotated <- data.frame(t(cpg_table))
str(cpg_table)
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              search = "random",
                              verboseIter = TRUE)
cpg_table_rotated$Age <- sample_table$Age
training <-train(Age ~ .,
      data = cpg_table_rotated,
      method = "glmnet",
      preProcess = c("center", "scale"),
      tuneLength = 25,
      trControl = train_control)
