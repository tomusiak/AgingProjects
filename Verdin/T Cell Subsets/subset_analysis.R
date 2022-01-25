#Grabs some useful scripts.
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")

#Sets location of data 
setwd("/home/atom/Desktop/Data/subset_data") #Sets directory.

#Libraries to import.
library(WGCNA)
library(missMethyl)
library(limma)
library(dplyr)
require("biomaRt")
library("RColorBrewer")
library(DMRcate)
library(readxl)
library(ggplot2)
library(GOfuncR)
library(readr)
library(plotrix)

clock_data <- read.csv("~/Desktop/Data/subset_data/clock_data.csv")
subset_metadata <- read.csv("~/Desktop/Data/subset_data/subset_metadata.csv")
all_data <- merge(clock_data,subset_metadata)
subset_data <- all_data[all_data$sabgal_sample==FALSE,]
subset_data$diff <- subset_data$age-subset_data$DNAmAge

summary <- summarySE(subset_data,"diff", "type")
summary$type <- c("Central Memory","Effector Memory","Naive","TEMRA")
summary$type <- factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=summary, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-18,18) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Age - Predicted Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age") 
