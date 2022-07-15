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
library(stringr)

sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
cpg_table <- data.table::fread("ClockConstruction/cpg_table.csv",header=TRUE)  %>% as.data.frame()
row.names(cpg_table) <- cpg_table$V1
cpg_table <- cpg_table[,-1]

# 
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
# 
# # Time for a monocyte data set.
# estupinan2022_unformatted_table <- read.table("Estupinan2022/GSE201752_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# estupinan2022_formatted_samples <- strsplit(estupinan2022_unformatted_table[1,2]," ")[[1]]
# 
# #Formatting ages.
# estupinan2022_unformatted_ages <- estupinan2022_unformatted_table[6,-1]
# estupinan2022_formatted_ages <- as.numeric(str_sub(estupinan2022_unformatted_ages, 6, 7))
# 
# #Formatting sex.
# estupinan2022_unformatted_sex <- estupinan2022_unformatted_table[5,-1]
# estupinan2022_formatted_sex <- str_sub(estupinan2022_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# estupinan2022_unformatted_condition <- estupinan2022_unformatted_table[4,-1]
# estupinan2022_formatted_condition <- str_sub(estupinan2022_unformatted_condition, 17, 20)
# estupinan2022_formatted_condition[estupinan2022_formatted_condition == "HD"] <- "Control"
# estupinan2022_formatted_misc <- estupinan2022_formatted_condition
# estupinan2022_formatted_condition[estupinan2022_formatted_condition != "Control"] <- "Giant Cell Arteritis"
# estupinan2022_IDs <- paste(rep("H",113),1:113,sep="")
# 
# estupinan2022_samples <- data.frame(ID=estupinan2022_formatted_samples,
#                                     Author=rep("Estupinan",113),
#                                     Year=rep(2022,113),
#                                     Tissue=rep("Blood",113),
#                                     CellType=rep("Monocytes",113),
#                                     Age=estupinan2022_formatted_ages,
#                                     Condition=estupinan2022_formatted_condition,
#                                     Sex=estupinan2022_formatted_sex,
#                                     DonorID=estupinan2022_IDs,
#                                     Misc=estupinan2022_formatted_misc)
# 
# estupinan2022_cpgs <- data.frame(read_table2("Estupinan2022/GSE201752_processed_data.txt", skip=4),row.names=1)
# colnames(estupinan2022_cpgs) <- estupinan2022_formatted_samples
# estupinan2022_cpgs <- estupinan2022_cpgs[rownames(estupinan2022_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(estupinan2022_cpgs),]
# 
# sample_table <- rbind(sample_table,estupinan2022_samples)
# cpg_table <- cbind(cpg_table,estupinan2022_cpgs)
# 
# 
# #Time for a blood data set.
# okereke2021_unformatted_table <- read.table("Okereke2021/GSE190540_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# okereke2021_unformatted_samples <- strsplit(okereke2021_unformatted_table[1,2]," ")[[1]]
# 
# #Formatting ages.
# okereke2021_unformatted_ages <- okereke2021_unformatted_table[4,-1]
# okereke2021_formatted_ages <- as.numeric(str_sub(okereke2021_unformatted_ages, 6, 7))
# 
# #Formatting sex.
# okereke2021_unformatted_sex <- okereke2021_unformatted_table[2,-1]
# okereke2021_formatted_sex <- str_sub(okereke2021_unformatted_sex, 9, 15)
# 
# #Formatting condition and miscellaneous information.
# okereke2021_unformatted_condition <- okereke2021_unformatted_table[3,-1]
# okereke2021_formatted_condition <- str_sub(okereke2021_unformatted_condition, 17, 25)
# okereke2021_formatted_condition[okereke2021_formatted_condition == "Case"] <- "Cognitive Impairment"
# 
# okereke2021_IDs <- paste(rep("I",90),1:90,sep="")
# 
# okereke2021_samples <- data.frame(ID=okereke2021_unformatted_samples,
#                                     Author=rep("Okereke",90),
#                                     Year=rep(2021,90),
#                                     Tissue=rep("Blood",90),
#                                     CellType=rep("PBMCs",90),
#                                     Age=okereke2021_formatted_ages,
#                                     Condition=okereke2021_formatted_condition,
#                                     Sex=okereke2021_formatted_sex,
#                                     DonorID=okereke2021_IDs,
#                                     Misc=rep(NA,90))
# 
# okereke2021_cpgs <- read.table("Okereke2021/GSE190540_series_matrix.txt",
#                                comment = "!",
#                                skip=5,
#                                fill=TRUE)
# okereke2021_cpgs <- data.frame(okereke2021_cpgs)
# okereke2021_cpgs <- okereke2021_cpgs[-c(1:5),]
# rownames(okereke2021_cpgs) <- okereke2021_cpgs$V1
# okereke2021_cpgs <- data.frame(okereke2021_cpgs[,-1])
# colnames(okereke2021_cpgs) <- okereke2021_unformatted_samples
# okereke2021_cpgs <- okereke2021_cpgs[rownames(okereke2021_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(okereke2021_cpgs),]
# 
# sample_table <- rbind(sample_table,okereke2021_samples)
# cpg_table <- cbind(cpg_table,okereke2021_cpgs)
# 
# # Skin data set.
# muse2021_unformatted_table <- read.table("Muse2021/GSE188593_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# muse2021_unformatted_samples <- muse2021_unformatted_table[1,-1]
# muse2021_formatted_samples <- t(muse2021_unformatted_samples)[,1]
# 
# #Formatting ages.
# muse2021_unformatted_ages <- muse2021_unformatted_table[3,-1]
# muse2021_formatted_ages <- as.numeric(str_sub(muse2021_unformatted_ages, 13, 19))
# 
# #Formatting sex.
# muse2021_unformatted_sex <- muse2021_unformatted_table[4,-1]
# muse2021_formatted_sex <- str_sub(muse2021_unformatted_sex, 6, 7)
# 
# muse2021_unformatted_donor <- muse2021_unformatted_table[2,-1]
# muse2021_unformatted_donor <- str_sub(muse2021_unformatted_donor, 10,15)
# muse2021_unformatted_donor <- as.numeric(factor(muse2021_unformatted_donor))
# muse2021_formatted_donor <- paste(rep("J",64),muse2021_unformatted_donor,sep="")
# 
# #Formatting condition and miscellaneous information.
# muse2021_unformatted_condition <- muse2021_unformatted_table[6,-1]
# muse2021_formatted_condition <- str_sub(muse2021_unformatted_condition, 17, 30)
# muse2021_formatted_condition[muse2021_formatted_condition == ""] <- "Control"
# 
# muse2021_IDs <- paste(rep("J",64),1:64,sep="")
# 
# muse2021_samples <- data.frame(ID=muse2021_formatted_samples,
#                                     Author=rep("Muse",64),
#                                     Year=rep(2021,64),
#                                     Tissue=rep("Skin",64),
#                                     CellType=rep("Epithelial",64),
#                                     Age=muse2021_formatted_ages,
#                                     Condition=muse2021_formatted_condition,
#                                     Sex=muse2021_formatted_sex,
#                                     DonorID=muse2021_IDs,
#                                     Misc=rep(NA,64))
# 
# muse2021_cpgs <- read.table("Muse2021/GSE188593_series_matrix.txt",
#                                comment = "!",
#                                skip=5,
#                                fill=TRUE)
# muse2021_cpgs <- data.frame(muse2021_cpgs)
# muse2021_cpgs <- muse2021_cpgs[-c(1:6),]
# rownames(muse2021_cpgs) <- muse2021_cpgs$V1
# muse2021_cpgs <- data.frame(muse2021_cpgs[,-1])
# colnames(muse2021_cpgs) <- muse2021_unformatted_samples
# muse2021_cpgs <- muse2021_cpgs[rownames(muse2021_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(muse2021_cpgs),]
# 
# sample_table <- rbind(sample_table,muse2021_samples)
# cpg_table <- cbind(cpg_table,muse2021_cpgs)

# #Blood data set.
# konigsberg2021_unformatted_table <- read.table("Konigsberg2021/GSE167202_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# konigsberg2021_formatted_samples <- strsplit(konigsberg2021_unformatted_table[1,2]," ")[[1]]
# 
# #Formatting ages.
# konigsberg2021_unformatted_ages <- konigsberg2021_unformatted_table[4,-1]
# konigsberg2021_formatted_ages <- as.numeric(str_sub(konigsberg2021_unformatted_ages, 5, 8))
# 
# #Formatting sex.
# konigsberg2021_unformatted_sex <- konigsberg2021_unformatted_table[3,-1]
# konigsberg2021_formatted_sex <- str_sub(konigsberg2021_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# konigsberg2021_unformatted_condition <- konigsberg2021_unformatted_table[2,-1]
# konigsberg2021_formatted_condition <- str_sub(konigsberg2021_unformatted_condition, 15, 30)
# konigsberg2021_formatted_condition[konigsberg2021_formatted_condition == "negative"] <- "Control"
# konigsberg2021_formatted_condition[konigsberg2021_formatted_condition == "positive"] <- "COVID"
# konigsberg2021_formatted_condition[konigsberg2021_formatted_condition == "other infection"] <-
#   "Respiratory Illness"
# 
# konigsberg2021_IDs <- paste(rep("K",525),1:525,sep="")
# 
# konigsberg2021_samples <- data.frame(ID=konigsberg2021_formatted_samples,
#                                     Author=rep("Konigsberg",525),
#                                     Year=rep(2021,525),
#                                     Tissue=rep("Blood",525),
#                                     CellType=rep("PBMCs",525),
#                                     Age=konigsberg2021_formatted_ages,
#                                     Condition=konigsberg2021_formatted_condition,
#                                     Sex=konigsberg2021_formatted_sex,
#                                     DonorID=konigsberg2021_IDs,
#                                     Misc=rep(NA,525))
# 
# konigsberg2021_cpgs <- data.table::fread("Konigsberg2021/GSE167202_ProcessedBetaValues.txt",
#                                          header=TRUE) %>%
#                       as.data.frame()
# konigsberg2021_cpgs <- konigsberg2021_cpgs[,-c(2:22)]
# konigsberg2021_sequence_ids <- konigsberg2021_unformatted_table[5,-1]
# konigsberg2021_sequence_ids <- unlist(konigsberg2021_sequence_ids)
# row.names(konigsberg2021_cpgs) <- konigsberg2021_cpgs$ID_REF
# konigsberg2021_cpgs <- konigsberg2021_cpgs[,-1]
# konigsberg2021_cpgs <- konigsberg2021_cpgs[,konigsberg2021_sequence_ids]
# colnames(konigsberg2021_cpgs) <- konigsberg2021_formatted_samples
# 
# konigsberg2021_cpgs <- konigsberg2021_cpgs[rownames(konigsberg2021_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(konigsberg2021_cpgs),]
# 
# sample_table <- rbind(sample_table,konigsberg2021_samples)
# cpg_table <- cbind(cpg_table,konigsberg2021_cpgs)

#Brain samples.
# haghighi2022_unformatted_table <- read.table("Haghighi2022/GSE191200_series_matrix.txt",
#                                                comment = "!",
#                                                skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# haghighi2022_formatted_samples <- strsplit(haghighi2022_unformatted_table[1,2]," ")[[1]]
# 
# #Formatting ages.
# haghighi2022_unformatted_ages <- haghighi2022_unformatted_table[4,-1]
# haghighi2022_formatted_ages <- as.numeric(str_sub(haghighi2022_unformatted_ages, 5, 8))
# 
# #Formatting sex.
# haghighi2022_unformatted_sex <- haghighi2022_unformatted_table[5,-1]
# haghighi2022_formatted_sex <- str_sub(haghighi2022_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# haghighi2022_unformatted_condition <- haghighi2022_unformatted_table[6,-1]
# haghighi2022_formatted_condition <- str_sub(haghighi2022_unformatted_condition, 12, 30)
# 
# haghighi2022_unformatted_donor <- haghighi2022_unformatted_table[2,-1]
# haghighi2022_unformatted_donor <- str_sub(haghighi2022_unformatted_donor,11,15)
# haghighi2022_unformatted_donor <- as.numeric(factor(haghighi2022_unformatted_donor,
#                                              labels=c(1:22)))
# haghighi2022_formatted_donor <- paste(rep("L",56),haghighi2022_unformatted_donor,sep="")
# 
# haghighi2022_samples <- data.frame(ID=haghighi2022_formatted_samples,
#                                      Author=rep("Haghighi",56),
#                                      Year=rep(2022,56),
#                                      Tissue=rep("Brain",56),
#                                      CellType=rep("Microglia",56),
#                                      Age=haghighi2022_formatted_ages,
#                                      Condition=haghighi2022_formatted_condition,
#                                      Sex=haghighi2022_formatted_sex,
#                                      DonorID=haghighi2022_formatted_donor,
#                                      Misc=rep(NA,56))
# 
# haghighi2022_cpgs <- read.table("Haghighi2022/GSE191200_series_matrix.txt",
#                                                                 comment = "!",
#                                                                 skip=7,
#                                                                 fill=TRUE)
# haghighi2022_cpgs <- haghighi2022_cpgs[-c(1:6),]
# row.names(haghighi2022_cpgs) <- haghighi2022_cpgs$V1
# haghighi2022_cpgs <- haghighi2022_cpgs[,-1]
# colnames(haghighi2022_cpgs) <- haghighi2022_cpgs[1,]
# haghighi2022_cpgs <- haghighi2022_cpgs[-1,]
# 
# haghighi2022_cpgs <- haghighi2022_cpgs[rownames(haghighi2022_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(haghighi2022_cpgs),]
# 
# sample_table <- rbind(sample_table,haghighi2022_samples)
# cpg_table <- cbind(cpg_table,haghighi2022_cpgs)

# #Colorectal samples.
# chen2021_unformatted_table <- read.table("Chen2021/GSE159898_series_matrix.txt",
#                                                comment = "!",
#                                                skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# chen2021_formatted_samples <- strsplit(chen2021_unformatted_table[1,2]," ")[[1]]
# 
# #Formatting ages.
# chen2021_unformatted_ages <- chen2021_unformatted_table[3,-1]
# chen2021_formatted_ages <- as.numeric(str_sub(chen2021_unformatted_ages, 5, 8))
# 
# #Formatting sex.
# chen2021_unformatted_sex <- chen2021_unformatted_table[2,-1]
# chen2021_formatted_sex <- str_sub(chen2021_unformatted_sex, 9, 15)
# 
# #Formatting condition and miscellaneous information.
# chen2021_unformatted_condition <- chen2021_unformatted_table[4,-1]
# chen2021_formatted_condition <- str_sub(chen2021_unformatted_condition, 14, 40)
# chen2021_formatted_condition[chen2021_formatted_condition == "normal colorectal tissue"] <-
#      "Control"
# chen2021_formatted_condition[chen2021_formatted_condition == "colorectal cancer tissue"] <-
#   "Colorectal Cancer"
# 
# chen2021_unformatted_donor <- chen2021_unformatted_table[5,-1]
# chen2021_unformatted_donor <- str_sub(chen2021_unformatted_donor,11,15)
# chen2021_unformatted_donor <- as.numeric(factor(chen2021_unformatted_donor,
#                                              labels=c(1:21)))
# chen2021_formatted_donor <- paste(rep("M",44),chen2021_unformatted_donor,sep="")
# 
# chen2021_samples <- data.frame(ID=chen2021_formatted_samples,
#                                      Author=rep("Chen",44),
#                                      Year=rep(2021,44),
#                                      Tissue=rep("Colorectal",44),
#                                      CellType=rep("Epithelial (Colorectal)",44),
#                                      Age=chen2021_formatted_ages,
#                                      Condition=chen2021_formatted_condition,
#                                      Sex=chen2021_formatted_sex,
#                                      DonorID=chen2021_formatted_donor,
#                                      Misc=rep(NA,44))
# 
# chen2021_cpgs <- read.table("Chen2021/GSE159898_series_matrix.txt",
#                                                                 comment = "!",
#                                                                 skip=7,
#                                                                 fill=TRUE)
# chen2021_cpgs <- chen2021_cpgs[-c(1:5),]
# row.names(chen2021_cpgs) <- chen2021_cpgs$V1
# chen2021_cpgs <- chen2021_cpgs[,-1]
# colnames(chen2021_cpgs) <- chen2021_cpgs[1,]
# chen2021_cpgs <- chen2021_cpgs[-1,]
# 
# chen2021_cpgs <- chen2021_cpgs[rownames(chen2021_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(chen2021_cpgs),]
# 
# sample_table <- rbind(sample_table,chen2021_samples)
# cpg_table <- cbind(cpg_table,chen2021_cpgs)

#Skeletal Muscle samples
# voisin2021_unformatted_table <- read.table("Voisin2021/GSE151407_series_matrix.txt",
#                                          comment = "!",
#                                          skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# voisin2021_formatted_samples <- strsplit(voisin2021_unformatted_table[1,2]," ")[[1]]
# 
# #Formatting ages.
# voisin2021_unformatted_ages <- voisin2021_unformatted_table[3,-1]
# voisin2021_formatted_ages <- as.numeric(str_sub(voisin2021_unformatted_ages, 5, 8))
# 
# #Formatting sex.
# voisin2021_unformatted_sex <- voisin2021_unformatted_table[4,-1]
# voisin2021_formatted_sex <- str_sub(voisin2021_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# voisin2021_unformatted_condition <- voisin2021_unformatted_table[2,-1]
# voisin2021_formatted_condition <- str_sub(voisin2021_unformatted_condition, 8, 10)
# voisin2021_formatted_condition[voisin2021_formatted_condition == "PRE"] <-
#   "Control"
# voisin2021_formatted_condition[voisin2021_formatted_condition != "Control"] <-
#   "HIIT"
# 
# voisin2021_unformatted_donor <- voisin2021_unformatted_table[2,-1]
# voisin2021_unformatted_donor <- str_sub(voisin2021_unformatted_donor,-9,-6)
# voisin2021_unformatted_donor <- sub("_", "", voisin2021_unformatted_donor)
# voisin2021_unformatted_donor <- as.numeric(factor(voisin2021_unformatted_donor,
#                                                 labels=c(1:25)))
# voisin2021_formatted_donor <- paste(rep("N",78),voisin2021_unformatted_donor,sep="")
# 
# voisin2021_samples <- data.frame(ID=voisin2021_formatted_samples,
#                                Author=rep("Voisin",78),
#                                Year=rep(2021,78),
#                                Tissue=rep("Skeletal Muscle",78),
#                                CellType=rep("Muscle Cells",78),
#                                Age=voisin2021_formatted_ages,
#                                Condition=voisin2021_formatted_condition,
#                                Sex=voisin2021_formatted_sex,
#                                DonorID=voisin2021_formatted_donor,
#                                Misc=rep(NA,78))
# 
# voisin2021_cpgs <- read.table("Voisin2021/GSE151407_series_matrix.txt",
#                             comment = "!",
#                             skip=7,
#                             fill=TRUE)
# voisin2021_cpgs <- voisin2021_cpgs[-c(1:4),]
# row.names(voisin2021_cpgs) <- voisin2021_cpgs$V1
# voisin2021_cpgs <- voisin2021_cpgs[,-1]
# colnames(voisin2021_cpgs) <- voisin2021_cpgs[1,]
# voisin2021_cpgs <- voisin2021_cpgs[-1,]
# 
# voisin2021_cpgs <- voisin2021_cpgs[rownames(voisin2021_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(voisin2021_cpgs),]
# 
# sample_table <- rbind(sample_table,voisin2021_samples)
# cpg_table <- cbind(cpg_table,voisin2021_cpgs)

# More brain. Note - I am somewhat apprehensive about this following dataset.

# fries2019_unformatted_table <- read.table("Fries2019/GSE129428_series_matrix.txt",
#                                          comment = "!",
#                                          skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# fries2019_formatted_samples <- strsplit(fries2019_unformatted_table[1,2]," ")[[1]]
# 
# #Formatting ages.
# fries2019_unformatted_ages <- fries2019_unformatted_table[4,-1]
# fries2019_formatted_ages <- as.numeric(str_sub(fries2019_unformatted_ages, 14, 19))
# 
# #Formatting sex.
# fries2019_unformatted_sex <- fries2019_unformatted_table[2,-1]
# fries2019_formatted_sex <- str_sub(fries2019_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# fries2019_unformatted_condition <- fries2019_unformatted_table[3,-1]
# fries2019_formatted_condition <- str_sub(fries2019_unformatted_condition, 8, 100)
# 
# fries2019_formatted_donor <- paste(rep("O",64),c(1:64),sep="")
# 
# fries2019_samples <- data.frame(ID=fries2019_formatted_samples,
#                                Author=rep("Fries",64),
#                                Year=rep(2019,64),
#                                Tissue=rep("Brain",64),
#                                CellType=rep("Brain Cells",64),
#                                Age=fries2019_formatted_ages,
#                                Condition=fries2019_formatted_condition,
#                                Sex=fries2019_formatted_sex,
#                                DonorID=fries2019_formatted_donor,
#                                Misc=rep(NA,64))
# 
# fries2019_cpgs <- read.table("Fries2019/GSE129428_series_matrix.txt",
#                             comment = "!",
#                             skip=7,
#                             fill=TRUE)
# fries2019_cpgs <- fries2019_cpgs[-c(1:4),]
# row.names(fries2019_cpgs) <- fries2019_cpgs$V1
# fries2019_cpgs <- fries2019_cpgs[,-1]
# colnames(fries2019_cpgs) <- fries2019_cpgs[1,]
# fries2019_cpgs <- fries2019_cpgs[-1,]
# 
# fries2019_cpgs <- fries2019_cpgs[rownames(fries2019_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(fries2019_cpgs),]
# 
# sample_table <- rbind(sample_table,fries2019_samples)
# cpg_table <- cbind(cpg_table,fries2019_cpgs)

#Samples from children.

# islam2018_unformatted_table <- read.table("Islam2018/GSE124366_series_matrix.txt",
#                                          comment = "!",
#                                          skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# islam2018_unformatted_samples <- islam2018_unformatted_table[6,]
# islam2018_formatted_samples <- str_sub(islam2018_unformatted_table[6,],1,100)[2:216]
# 
# #Formatting ages.
# islam2018_unformatted_ages <- islam2018_unformatted_table[4,-1]
# islam2018_formatted_ages <- as.numeric(str_sub(islam2018_unformatted_ages, 35,100))
# 
# #Formatting sex.
# islam2018_unformatted_sex <- islam2018_unformatted_table[3,-1]
# islam2018_formatted_sex <- str_sub(islam2018_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# islam2018_formatted_condition <- rep("Control",215)
# 
# islam2018_unformatted_donor <- islam2018_unformatted_table[1,-1]
# islam2018_unformatted_donor <- str_sub(islam2018_unformatted_donor,1,8)
# islam2018_unformatted_donor <- as.numeric(factor(islam2018_unformatted_donor,
#                                                 labels=c(1:201)))
# islam2018_formatted_donor <- paste(rep("P",201),islam2018_unformatted_donor,sep="")
# 
# islam2018_unformatted_celltype <- islam2018_unformatted_table[2,-1]
# islam2018_formatted_celltype <- str_sub(islam2018_unformatted_celltype, 9, 20)
# 
# islam2018_unformatted_tissue <- islam2018_formatted_celltype
# islam2018_unformatted_tissue[islam2018_unformatted_tissue == "PBMC"] <- "Blood"
# islam2018_formatted_tissue <- islam2018_unformatted_tissue
# 
# islam2018_samples <- data.frame(ID=islam2018_formatted_samples,
#                                Author=rep("Islam",215),
#                                Year=rep(2018,215),
#                                Tissue=islam2018_formatted_tissue,
#                                CellType=islam2018_formatted_celltype,
#                                Age=islam2018_formatted_ages,
#                                Condition=islam2018_formatted_condition,
#                                Sex=islam2018_formatted_sex,
#                                DonorID=islam2018_formatted_donor,
#                                Misc=rep(NA,215))
# 
# islam2018_cpgs <- read.table("Islam2018/GSE124366_series_matrix.txt",
#                             comment = "!",
#                             skip=6,
#                             fill=TRUE)
# islam2018_cpgs <- islam2018_cpgs[-c(1:6),]
# row.names(islam2018_cpgs) <- islam2018_cpgs$V1
# islam2018_cpgs <- islam2018_cpgs[,-1]
# colnames(islam2018_cpgs) <- islam2018_formatted_samples
# 
# islam2018_cpgs <- islam2018_cpgs[rownames(islam2018_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(islam2018_cpgs),]
# 
# sample_table <- rbind(sample_table,islam2018_samples)
# cpg_table <- cbind(cpg_table,islam2018_cpgs)
# 

#Samples from breast tissue.

# johnson2016_unformatted_table <- read.table("Johnson2016/GSE88883_series_matrix.txt",
#                                           comment = "!",
#                                           skip=0,fill=TRUE,nrows = 4)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# johnson2016_unformatted_samples <- johnson2016_unformatted_table[4,]
# johnson2016_formatted_samples <- str_sub(johnson2016_unformatted_samples,1,100)[2:101]
# 
# #Formatting ages.
# johnson2016_unformatted_ages <- johnson2016_unformatted_table[2,-1]
# johnson2016_formatted_ages <- as.numeric(str_sub(johnson2016_unformatted_ages, 14,100))
# 
# #Formatting sex.
# johnson2016_unformatted_sex <- johnson2016_unformatted_table[1,-1]
# johnson2016_formatted_sex <- str_sub(johnson2016_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# johnson2016_formatted_condition <- rep("Control",100)
# 
# johnson2016_formatted_donor <- paste(rep("Q",100),c(1:100),sep="")
# 
# johnson2016_samples <- data.frame(ID=johnson2016_formatted_samples,
#                                 Author=rep("Johnson",100),
#                                 Year=rep(2016,100),
#                                 Tissue=rep("Breast",100),
#                                 CellType=rep("Breast",100),
#                                 Age=johnson2016_formatted_ages,
#                                 Condition=johnson2016_formatted_condition,
#                                 Sex=johnson2016_formatted_sex,
#                                 DonorID=johnson2016_formatted_donor,
#                                 Misc=rep(NA,100))
# 
# johnson2016_cpgs <- read.table("Johnson2016/GSE88883_series_matrix.txt",
#                              comment = "!",
#                              skip=4,
#                              fill=TRUE)
# johnson2016_cpgs <- johnson2016_cpgs[-c(1:4),]
# row.names(johnson2016_cpgs) <- johnson2016_cpgs$V1
# johnson2016_cpgs <- johnson2016_cpgs[,-1]
# colnames(johnson2016_cpgs) <- johnson2016_formatted_samples
# 
# johnson2016_cpgs <- johnson2016_cpgs[rownames(johnson2016_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(johnson2016_cpgs),]
# 
# sample_table <- rbind(sample_table,johnson2016_samples)
# cpg_table <- cbind(cpg_table,johnson2016_cpgs)

#Samples from liver.
# 
# horvath2014_unformatted_table <- read.table("Horvath2014/GSE61258_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 4)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# horvath2014_unformatted_samples <- horvath2014_unformatted_table[4,]
# horvath2014_formatted_samples <- str_sub(horvath2014_unformatted_samples,1,100)[2:80]
# 
# #Formatting ages.
# horvath2014_unformatted_ages <- horvath2014_unformatted_table[2,-1]
# horvath2014_formatted_ages <- as.numeric(str_sub(horvath2014_unformatted_ages, 5,100))
# 
# #Formatting sex.
# horvath2014_unformatted_sex <- horvath2014_unformatted_table[1,-1]
# horvath2014_formatted_sex <- str_sub(horvath2014_unformatted_sex, 6, 15)
# 
# #Formatting condition and miscellaneous information.
# horvath2014_formatted_condition <- rep("Control",79)
# 
# horvath2014_formatted_donor <- paste(rep("R",79),c(1:79),sep="")
# 
# horvath2014_samples <- data.frame(ID=horvath2014_formatted_samples,
#                                   Author=rep("Horvath",79),
#                                   Year=rep(2014,79),
#                                   Tissue=rep("Liver",79),
#                                   CellType=rep("Liver",79),
#                                   Age=horvath2014_formatted_ages,
#                                   Condition=horvath2014_formatted_condition,
#                                   Sex=horvath2014_formatted_sex,
#                                   DonorID=horvath2014_formatted_donor,
#                                   Misc=rep(NA,79))
# 
# horvath2014_cpgs <- read.table("Horvath2014/GSE61258_series_matrix.txt",
#                                comment = "!",
#                                skip=4,
#                                fill=TRUE)
# horvath2014_cpgs <- horvath2014_cpgs[-c(1:4),]
# row.names(horvath2014_cpgs) <- horvath2014_cpgs$V1
# horvath2014_cpgs <- horvath2014_cpgs[,-1]
# colnames(horvath2014_cpgs) <- horvath2014_formatted_samples
# 
# horvath2014_cpgs <- horvath2014_cpgs[rownames(horvath2014_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(horvath2014_cpgs),]
# 
# sample_table <- rbind(sample_table,horvath2014_samples)
# cpg_table <- cbind(cpg_table,horvath2014_cpgs)

#Samples from sperm tissue.

# pilsner2022_unformatted_table <- read.table("Pilsner2022/GSE185445_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 5)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# pilsner2022_unformatted_samples <- pilsner2022_unformatted_table[1,-1]
# pilsner2022_formatted_samples <- str_sub(pilsner2022_unformatted_samples,1,100)[1:379]
# 
# #Formatting ages.
# pilsner2022_unformatted_ages <- pilsner2022_unformatted_table[3,-1]
# pilsner2022_formatted_ages <- as.numeric(str_sub(pilsner2022_unformatted_ages, 8,100))
# 
# #Formatting sex.
# pilsner2022_formatted_sex <- rep("Male",379)
# 
# #Formatting condition and miscellaneous information.
# pilsner2022_formatted_condition <- rep("Control",379)
# 
# pilsner2022_formatted_donor <- paste(rep("T",379),c(1:379),sep="")
# 
# pilsner2022_samples <- data.frame(ID=pilsner2022_formatted_samples,
#                                   Author=rep("Pilsner",379),
#                                   Year=rep(2022,379),
#                                   Tissue=rep("Semen",379),
#                                   CellType=rep("Sperm",379),
#                                   Age=pilsner2022_formatted_ages,
#                                   Condition=pilsner2022_formatted_condition,
#                                   Sex=pilsner2022_formatted_sex,
#                                   DonorID=pilsner2022_formatted_donor,
#                                   Misc=rep(NA,379))
# 
# pilsner2022_unformatted_basenames <- pilsner2022_unformatted_table[2,-1]
# pilsner2022_formatted_basenames <- str_sub(pilsner2022_unformatted_basenames,11,100)
# 
# pilsner2022_cpgs <- data.table::fread("Pilsner2022/GSE185445_processed_methylation_data_matrix.txt",
#                                       header=FALSE)  %>% 
#                                       as.data.frame()
# 
# pilsner2022_colnames <- data.frame(read_table2("Pilsner2022/colnames.csv", 
#                         col_names = FALSE))
# pilsner2022_colnames <- substr(pilsner2022_colnames,2,20)
# 
# row.names(pilsner2022_cpgs) <- pilsner2022_cpgs[,1]
# pilsner2022_cpgs <- pilsner2022_cpgs[,-1]
# colnames(pilsner2022_cpgs) <- pilsner2022_colnames
# pilsner2022_cpgs <- pilsner2022_cpgs[,pilsner2022_formatted_basenames]
# 
# pilsner2022_cpgs <- pilsner2022_cpgs[rownames(pilsner2022_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(pilsner2022_cpgs),]
# colnames(pilsner2022_cpgs) <- pilsner2022_formatted_samples
# 
# sample_table <- rbind(sample_table,pilsner2022_samples)
# cpg_table <- cbind(cpg_table,pilsner2022_cpgs)

#Samples from nasal epithelia

# davalos2022_unformatted_table <- read.table("Davalos2022/GSE193879_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 4)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# davalos2022_unformatted_samples <- davalos2022_unformatted_table[1,]
# davalos2022_formatted_samples <- str_sub(davalos2022_unformatted_samples,1,62)[2:128]
# 
# #Formatting ages.
# davalos2022_unformatted_ages <- davalos2022_unformatted_table[4,-1]
# davalos2022_formatted_ages <- as.numeric(str_sub(davalos2022_unformatted_ages, 13,100))
# 
# #Formatting sex.
# davalos2022_unformatted_sex <- davalos2022_unformatted_table[3,-1]
# davalos2022_formatted_sex <- str_sub(davalos2022_unformatted_sex, 17, 30)
# 
# #Formatting condition and miscellaneous information.
# davalos2022_unformatted_condition <- davalos2022_unformatted_table[2,-1]
# davalos2022_formatted_condition <- str_sub(davalos2022_unformatted_condition, 8,39)
# 
# davalos2022_formatted_donor <- paste(rep("U",127),c(1:127),sep="")
# 
# davalos2022_samples <- data.frame(ID=davalos2022_formatted_samples,
#                                   Author=rep("Davalos",127),
#                                   Year=rep(2022,127),
#                                   Tissue=rep("Blood",127),
#                                   CellType=rep("Blood",127),
#                                   Age=davalos2022_formatted_ages,
#                                   Condition=davalos2022_formatted_condition,
#                                   Sex=davalos2022_formatted_sex,
#                                   DonorID=davalos2022_formatted_donor,
#                                   Misc=rep(NA,127))
# 
# davalos2022_samples <- davalos2022_samples[!is.na(davalos2022_samples$Age),]
# 
# davalos2022_cpgs <- data.table::fread("Davalos2022/GSE193879_Matrix_processed.csv",
#                                                       header=FALSE)  %>% 
#                                                       as.data.frame()
# row.names(davalos2022_cpgs) <- davalos2022_cpgs$V1
# davalos2022_cpgs <- davalos2022_cpgs[-1,-1]
# davalos2022_cpgs <- davalos2022_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),127)]
# colnames(davalos2022_cpgs) <- davalos2022_formatted_samples
# davalos2022_cpgs <- davalos2022_cpgs[,colnames(davalos2022_cpgs) %in% davalos2022_samples$ID]
# 
# davalos2022_cpgs <- davalos2022_cpgs[rownames(davalos2022_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(davalos2022_cpgs),]
# 
# sample_table <- rbind(sample_table,davalos2022_samples)
# cpg_table <- cbind(cpg_table,davalos2022_cpgs)


# hannon2021_unformatted_table <- read.table("Hannon2021/GSE152026_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 6)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# hannon2021_unformatted_samples <- hannon2021_unformatted_table[2,]
# hannon2021_formatted_samples <- str_sub(hannon2021_unformatted_samples,1,62)[2:935]
# 
# #Formatting ages.
# hannon2021_unformatted_ages <- hannon2021_unformatted_table[5,-1]
# hannon2021_formatted_ages <- as.numeric(str_sub(hannon2021_unformatted_ages, 5,100))
# 
# #Formatting sex.
# hannon2021_unformatted_sex <- hannon2021_unformatted_table[4,-1]
# hannon2021_formatted_sex <- str_sub(hannon2021_unformatted_sex, 6, 30)
# 
# #Formatting condition and miscellaneous information.
# hannon2021_unformatted_condition <- hannon2021_unformatted_table[3,-1]
# hannon2021_formatted_condition <- str_sub(hannon2021_unformatted_condition, 12,39)
# hannon2021_formatted_condition[hannon2021_formatted_condition=="Case"] <- "Schizophrenia"
# 
# hannon2021_formatted_donor <- paste(rep("V",934),c(1:934),sep="")
# 
# hannon2021_samples <- data.frame(ID=hannon2021_formatted_samples,
#                                   Author=rep("Hannon",934),
#                                   Year=rep(2021,934),
#                                   Tissue=rep("Blood",934),
#                                   CellType=rep("Blood",934),
#                                   Age=hannon2021_formatted_ages,
#                                   Condition=hannon2021_formatted_condition,
#                                   Sex=hannon2021_formatted_sex,
#                                   DonorID=hannon2021_formatted_donor,
#                                   Misc=rep(NA,934))
# 
# hannon2021_samples <- hannon2021_samples[!is.na(hannon2021_samples$Age),]
# hannon2021_sample_names <- str_sub(hannon2021_unformatted_table[1,-1],1,19)
# 
# hannon2021_cpgs <- data.table::fread("Hannon2021/GSE152026_EUGEI_processed_signals.csv",
#                                       header=FALSE)  %>% 
#   as.data.frame()
# row.names(hannon2021_cpgs) <- hannon2021_cpgs$V1
# hannon2021_cpgs <- hannon2021_cpgs[-1,-1]
# hannon2021_cpgs <- hannon2021_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),934)]
# colnames(hannon2021_cpgs) <- hannon2021_formatted_samples
# hannon2021_cpgs <- hannon2021_cpgs[,colnames(hannon2021_cpgs) %in% hannon2021_samples$ID]
# hannon2021_names <- read_csv("Hannon2021/target.txt", 
#                              col_names = FALSE)
# hannon2021_names <- str_sub(unlist(hannon2021_names),1,500)
# hannon2021_names <- hannon2021_names[2:1869]
# hannon2021_names <- hannon2021_names[rep(c(rep(TRUE, 2- 1), FALSE),934)]
# colnames(hannon2021_cpgs) <- hannon2021_names
# 
# hannon2021_cpgs <- hannon2021_cpgs[rownames(hannon2021_cpgs) %in% rownames(cpg_table),]
# 
# hannon2021_changed_cpgs <- hannon2021_cpgs[,hannon2021_sample_names]
# colnames(hannon2021_changed_cpgs) <- hannon2021_formatted_samples
# hannon2021_changed_cpgs <- hannon2021_changed_cpgs[,colnames(hannon2021_changed_cpgs) %in% 
#                                                      hannon2021_samples$ID]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(hannon2021_changed_cpgs),]
# 
# sample_table <- rbind(sample_table,hannon2021_samples)
# cpg_table <- cbind(cpg_table,hannon2021_changed_cpgs)

# martino2018_unformatted_table <- read.table("Martino2018/GSE114135-GPL23976_series_matrix.txt",
#                                             comment = "!",
#                                             skip=0,fill=TRUE,nrows = 5)
# 
# #Formatting sample names, as the header text file is not arranged in a particularly easy way.
# martino2018_unformatted_samples <- martino2018_unformatted_table[5,-1]
# martino2018_formatted_samples <- str_sub(martino2018_unformatted_samples,1,62)[1:205]
# 
# #Formatting ages.
# martino2018_unformatted_ages <- martino2018_unformatted_table[2,-1]
# martino2018_formatted_ages <- as.numeric(str_sub(martino2018_unformatted_ages, 5,100))
# 
# #Formatting sex.
# martino2018_unformatted_sex <- martino2018_unformatted_table[3,-1]
# martino2018_formatted_sex <- str_sub(martino2018_unformatted_sex, 6, 30)
# 
# #Formatting condition and miscellaneous information.
# martino2018_unformatted_condition <- martino2018_unformatted_table[4,-1]
# martino2018_formatted_condition <- str_sub(martino2018_unformatted_condition, -1,-1)
# martino2018_formatted_condition[martino2018_formatted_condition=="1"] <- "Stimulated"
# martino2018_formatted_condition[martino2018_formatted_condition=="0"] <- "Control"
# 
# martino2018_unformatted_donor <- martino2018_unformatted_table[1,-1]
# martino2018_unformatted_donor <- str_sub(martino2018_unformatted_donor,1,6)
# martino2018_unformatted_donor <- as.numeric(factor(martino2018_unformatted_donor))
# martino2018_formatted_donor <- paste(rep("W",205),martino2018_unformatted_donor,sep="")
# 
# martino2018_samples <- data.frame(ID=martino2018_formatted_samples,
#                                   Author=rep("Martino",205),
#                                   Year=rep(2018,205),
#                                   Tissue=rep("Blood",205),
#                                   CellType=rep("CD4+",205),
#                                   Age=martino2018_formatted_ages,
#                                   Condition=martino2018_formatted_condition,
#                                   Sex=martino2018_formatted_sex,
#                                   DonorID=martino2018_formatted_donor,
#                                   Misc=rep(NA,205))
# 
# martino2018_cpgs <- read.table("Martino2018/GSE114135-GPL23976_series_matrix.txt",
#                               comment = "!",
#                               skip=5,
#                               fill=TRUE)
# 
# martino2018_cpgs <- martino2018_cpgs[-c(1:5),]
# rownames(martino2018_cpgs) <- str_sub(martino2018_cpgs[,1],1,100)
# martino2018_cpgs <- martino2018_cpgs[,-1]
# colnames(martino2018_cpgs) <- martino2018_formatted_samples
# 
# martino2018_cpgs <- martino2018_cpgs[rownames(martino2018_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(martino2018_cpgs),]
# 
# sample_table <- rbind(sample_table,martino2018_samples)
# cpg_table <- cbind(cpg_table,martino2018_cpgs)


thompson2020_unformatted_table <- read.table("Thompson2020/GSE146376_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 7)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
thompson2020_unformatted_samples <- thompson2020_unformatted_table[7,-1]
thompson2020_formatted_samples <- str_sub(thompson2020_unformatted_samples,1,62)[1:280]

#Formatting ages.
thompson2020_unformatted_ages <- thompson2020_unformatted_table[3,-1]
thompson2020_formatted_ages <- as.numeric(str_sub(thompson2020_unformatted_ages, 5,100))

#Formatting sex.
thompson2020_unformatted_sex <- thompson2020_unformatted_table[4,-1]
thompson2020_formatted_sex <- str_sub(thompson2020_unformatted_sex, 6, 30)

#Formatting condition and miscellaneous information.
thompson2020_unformatted_condition <- thompson2020_unformatted_table[2,-1]
thompson2020_formatted_condition <- str_sub(thompson2020_unformatted_condition, 12, 30)
thompson2020_formatted_condition[thompson2020_formatted_condition=="Vehicle"] <- "Control"

thompson2020_unformatted_donor <- thompson2020_unformatted_table[1,-1]
thompson2020_unformatted_donor <- str_sub(thompson2020_unformatted_donor,8,10)
thompson2020_formatted_donor <- paste(rep("X",280),thompson2020_unformatted_donor,sep="")

thompson2020_unformatted_arrangement <- thompson2020_unformatted_table[5,-1]
thompson2020_formatted_arrangement <- str_sub(thompson2020_unformatted_arrangement,1,100)

thompson2020_samples <- data.frame(ID=thompson2020_formatted_samples,
                                  Author=rep("Thompson",280),
                                  Year=rep(2020,280),
                                  Tissue=rep("Lung",56),
                                  CellType=rep("Airway Smooth Muscle",56),
                                  Age=thompson2020_formatted_ages,
                                  Condition=thompson2020_formatted_condition,
                                  Sex=thompson2020_formatted_sex,
                                  DonorID=thompson2020_formatted_donor,
                                  Misc=rep(NA,280))

thompson2020_cpgs <- data.table::fread("Thompson2020/GSE146376_ProcessedMethData_70_ForGEO.csv",
                                      header=FALSE)  %>%
                                      as.data.frame()
row.names(thompson2020_cpgs) <- thompson2020_cpgs$V1
thompson2020_cpgs <- thompson2020_cpgs[,-1]
thompson2020_cpgs <- thompson2020_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),280)]

thompson2020_colnames <- read_csv("Thompson2020/colnames.csv", 
                     col_names = FALSE)
thompson2020_colnames <- str_sub(unlist(thompson2020_colnames[1,]),1,100)
thompson2020_colnames <- thompson2020_colnames[-1]
thompson2020_colnames <- thompson2020_colnames[rep(c(rep(TRUE, 2- 1), FALSE),280)]
colnames(thompson2020_cpgs) <- thompson2020_colnames

colnames(thompson2020_cpgs) <- thompson2020_formatted_samples

thompson2020_cpgs <- thompson2020_cpgs[rownames(thompson2020_cpgs) %in% rownames(cpg_table),]
cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(thompson2020_cpgs),]

sample_table <- rbind(sample_table,thompson2020_samples)
cpg_table <- cbind(cpg_table,thompson2020_cpgs)

###########################################################################

cpg_table <- na.omit(cpg_table)

healthy_samples <- sample_table[sample_table$Condition=="Control",1]
healthy_sample_table <- sample_table[sample_table$Condition=="Control",]
healthy_cpg_table <- cpg_table[,colnames(cpg_table) %in% healthy_samples]

data.table::fwrite(cpg_table,"ClockConstruction/cpg_table.csv", 
                   row.names = TRUE)
data.table::fwrite(sample_table,"ClockConstruction/sample_table.csv",
                   row.names = TRUE)
data.table::fwrite(healthy_cpg_table,"ClockConstruction/healthy_cpg_table.csv", 
                   row.names = TRUE)
data.table::fwrite(healthy_sample_table,"ClockConstruction/healthy_sample_table.csv",
                   row.names = TRUE)
