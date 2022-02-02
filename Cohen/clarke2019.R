#keep_mitogenes <- c("RNR1","RNR2", "ND1","ND2","CO1","CO2","ATP8",
#                    "ATP6","CO3","ND3","ND4L","ND4","ND5","ND6","CYB")

library("BiocManager")
library("S4Vectors")
library("Biostrings")
library("SummarizedExperiment")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("Biobase")
library("AnnotationDbi")
library("GenomicFeatures")
library(dplyr)
library(ggplot2)
library("RColorBrewer")
library(Rsubread)
library(ggrepel)

#Loading in auxiliary code.
setwd("/home/atom/Desktop/AgingProjects/Cohen/")
source("mdpseq_background_correction.R")

setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")
setwd("/home/atom/Desktop/Data/clarke2019")
bam_list <- makeBAMS(".",FALSE)
mitogene_counts <- getCountsMitochondrial(bam_list,FALSE)
mdp_counts <- getCountsMDP(bam_list,FALSE)
write.csv(mitogene_counts,"mitogene_counts.csv")
write.csv(mdp_counts,"mdp_counts.csv")

mdp_counts <- read.csv("mdp_counts.csv", row.names=1)
mitogene_counts <- read.csv("mitogene_counts.csv",row.names=1)
sun_samplesheet <- read.csv("sun_samplesheet.csv")
colnames(mdp_counts) <- sun_samplesheet$sample
colnames(mitogene_counts) <- sun_samplesheet$sample
sun_samplesheet$diagnosis <- as.factor(sun_samplesheet$diagnosis)