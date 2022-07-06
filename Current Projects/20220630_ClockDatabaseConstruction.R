#In this program I will construct a database of DNA methylation data, consisting of 450K and EPIC
# chip results. I will construct this database as having two segments - "cpg_table" consisting
# of cpgs in rows and samples in columns, and "sample_table" consisting of samples in rows and
# metadata in columns. The possible cpgs will be filtered on those that do not change with 
# T cell differentiation.

#NOTE: EPIC samples are 450k samples and vice versa.

source("AgingProjects/Useful Scripts/generally_useful.R") #Helper functions

#Packages
setwd("Data/") #Sets directory.
library(readr)
library(tidyr)
library(dplyr)
library(methylumi)
library(minfi)

# #Reading in the 450K dataset from Magnaye 2022 et al.
# magnaye2022450k_unformatted_table <- read.table("Magnaye2022/GSE201872-GPL21145_series_matrix.txt",
#                                             comment = "!",
#                                             skip=5,fill=TRUE,nrows = 5)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# magnaye2022450k_unformatted_samples <- magnaye2022450k_unformatted_table[2,-1]
# magnaye2022450k_formatted_samples <- t(magnaye2022450k_unformatted_samples)[,1]
# 
# #Formatting ages.
# magnaye2022450k_unformatted_ages <- magnaye2022450k_unformatted_table[3,-1]
# magnaye2022450k_formatted_ages <- as.numeric(str_sub(magnaye2022450k_unformatted_ages, 6, 7))
# 
# #Formatting sex.
# magnaye2022450k_unformatted_sex <- magnaye2022450k_unformatted_table[4,-1]
# magnaye2022450k_formatted_sex <- str_sub(magnaye2022450k_unformatted_sex, 6, 15)
# 
# #Formatting condition.
# magnaye2022450k_unformatted_condition <- magnaye2022450k_unformatted_table[5,-1]
# magnaye2022450k_formatted_condition <- str_sub(magnaye2022450k_unformatted_condition, 9, 18)
# 
# #Looks like the values provided are heavily pre-processed m-values. To ensure consistency,
# # I will process raw data to generate normalized beta values.
# magnaye2022450k_list_of_files<-list.files(file.path("Magnaye2022/raw_data"))
# magnaye2022450k_list_of_files <- magnaye2022450k_list_of_files[231:length(magnaye2022450k_list_of_files)]
# magnaye2022450k_parsed_list_of_files <- substr(magnaye2022450k_list_of_files,1,30)
# magnaye2022450k_parsed_list_of_files <- unique(magnaye2022450k_parsed_list_of_files)
# magnaye2022450k_samplesheet <- data.frame(Sample=magnaye2022450k_formatted_samples,each=2,
#                                        Ages = magnaye2022450k_formatted_ages,each=2,
#                                        Sex= magnaye2022450k_formatted_sex,each=2,
#                                        Condition = magnaye2022450k_formatted_condition,each=2,
#                                        Basename = magnaye2022450k_parsed_list_of_files)
# setwd("Magnaye2022/raw_data")
# magnaye2022450k_RGSet <- read.metharray.exp(targets = magnaye2022450k_samplesheet)
# magnaye2022450k_MSet <- preprocessSWAN(magnaye2022450k_RGSet)
# magnaye2022450k_cpgs <- getBeta(magnaye2022450k_MSet)
# setwd("..")
# setwd("..")
# colnames(magnaye2022450k_cpgs) <- magnaye2022450k_formatted_samples
# magnaye2022450k_cpgs <- data.frame(magnaye2022450k_cpgs)
# 
# #Reading in the EPIC dataset from Magnaye 2022 et al.
# magnaye2022EPIC_unformatted_table <- read.table("Magnaye2022/GSE201872-GPL13534_series_matrix.txt",
#                                                 comment = "!",
#                                                 skip=5,fill=TRUE,nrows = 5)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# magnaye2022EPIC_unformatted_samples <- magnaye2022EPIC_unformatted_table[1,-1]
# magnaye2022EPIC_formatted_samples <- t(magnaye2022EPIC_unformatted_samples)[,1]
# 
# #Formatting ages.
# magnaye2022EPIC_unformatted_ages <- magnaye2022EPIC_unformatted_table[2,-1]
# magnaye2022EPIC_formatted_ages <- as.numeric(str_sub(magnaye2022EPIC_unformatted_ages, 6, 7))
# 
# #Formatting sex.
# magnaye2022EPIC_unformatted_sex <- magnaye2022EPIC_unformatted_table[3,-1]
# magnaye2022EPIC_formatted_sex <- str_sub(magnaye2022EPIC_unformatted_sex, 6, 15)
# 
# #Formatting condition.
# magnaye2022EPIC_unformatted_condition <- magnaye2022EPIC_unformatted_table[4,-1]
# magnaye2022EPIC_formatted_condition <- str_sub(magnaye2022EPIC_unformatted_condition, 9, 18)
# 
# #Looks like the values provided are heavily pre-processed m-values. To ensure consistency,
# # I will process raw data to generate normalized beta values.
# magnaye2022EPIC_list_of_files<-list.files(file.path("Magnaye2022/raw_data"))
# magnaye2022EPIC_list_of_files <- magnaye2022EPIC_list_of_files[7:230]
# magnaye2022EPIC_parsed_list_of_files <- substr(magnaye2022EPIC_list_of_files,1,28)
# magnaye2022EPIC_parsed_list_of_files <- unique(magnaye2022EPIC_parsed_list_of_files)
# magnaye2022EPIC_samplesheet <- data.frame(Sample=magnaye2022EPIC_formatted_samples,each=2,
#                                           Ages = magnaye2022EPIC_formatted_ages,each=2,
#                                           Sex= magnaye2022EPIC_formatted_sex,each=2,
#                                           Condition = magnaye2022EPIC_formatted_condition,each=2,
#                                           Basename = magnaye2022EPIC_parsed_list_of_files)
# setwd("Magnaye2022/raw_data")
# magnaye2022EPIC_RGSet <- read.metharray.exp(targets = magnaye2022EPIC_samplesheet)
# magnaye2022EPIC_MSet <- preprocessSWAN(magnaye2022EPIC_RGSet)
# magnaye2022EPIC_cpgs <- getBeta(magnaye2022EPIC_MSet)
# setwd("..")
# setwd("..")
# colnames(magnaye2022EPIC_cpgs) <- magnaye2022EPIC_formatted_samples
# magnaye2022EPIC_cpgs <- data.frame(magnaye2022EPIC_cpgs)
# magnaye2022EPIC_cpgs <- magnaye2022EPIC_cpgs[rownames(magnaye2022EPIC_cpgs) %in%
#                                                rownames(magnaye2022450k_cpgs),]
# magnaye2022450k_cpgs <- magnaye2022450k_cpgs[rownames(magnaye2022450k_cpgs) %in%
#                                                rownames(magnaye2022EPIC_cpgs),]
# 
# #Initialize sample table.
# sample_table <- data.frame(ID=character(),
#                            Author=character(),
#                            Year=integer(),
#                            Tissue=character(),
#                            CellType=character(),
#                            Age=integer(),
#                            Condition=character(),
#                            Sex=character(),
#                            DonorID=character(),
#                            Misc=character())
# magnayefirst_IDs <- paste(rep("F",30),1:30,sep="")
# magnayefirst_samples <- data.frame(ID=magnaye2022450k_formatted_samples,
#                                    Author=rep("Magnaye",30),
#                                    Year=rep(2022,30),
#                                    Tissue=rep("Bronchi",30),
#                                    CellType=rep("Epithelial",30),
#                                    Age=magnaye2022450k_formatted_ages,
#                                    Condition=magnaye2022450k_formatted_condition,
#                                    Sex=magnaye2022450k_formatted_sex,
#                                    DonorID=magnayefirst_IDs,
#                                    Misc=rep("",30))
# magnayesecond_IDs <- paste(rep("G",112),1:112,sep="")
# magnayesecond_samples <- data.frame(ID=magnaye2022EPIC_formatted_samples,
#                                    Author=rep("Magnaye",112),
#                                    Year=rep(2022,112),
#                                    Tissue=rep("Bronchi",112),
#                                    CellType=rep("Epithelial",112),
#                                    Age=magnaye2022EPIC_formatted_ages,
#                                    Condition=magnaye2022EPIC_formatted_condition,
#                                    Sex=magnaye2022EPIC_formatted_sex,
#                                    DonorID=magnayesecond_IDs,
#                                    Misc=rep("",112))
# sample_table <- rbind(magnayefirst_samples,magnayesecond_samples)
# cpg_table <- cbind(magnaye2022450k_cpgs,magnaye2022EPIC_cpgs)
# 
# write.csv(sample_table,"ClockConstruction/sample_table.csv")
# write.csv(cpg_table,"ClockConstruction/cpg_table.csv")

sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
cpg_table <- read.csv("ClockConstruction/cpg_table.csv",row.names=1)