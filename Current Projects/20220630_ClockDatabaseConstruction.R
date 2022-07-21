#In this program I will construct a database of DNA methylation data, consisting of 450K and EPIC
# chip results. I will construct this database as having two segments - "cpg_table" consisting
# of cpgs in rows and samples in columns, and "sample_table" consisting of samples in rows and
# metadata in columns. The possible cpgs will be filtered on those that do not change with 
# T cell differentiation.

source("AgingProjects/Useful Scripts/generally_useful.R") #Helper functions

#Packages
setwd("Data/") #Sets directory.
library(readr)
library(tidyr)
library(dplyr)
library(methylumi)
library(minfi)
library(stringr)
library(data.table)

sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
cpg_table <- data.table::fread("ClockConstruction/cpg_table.csv",header=TRUE)  %>% as.data.frame()
row.names(cpg_table) <- cpg_table$V1
cpg_table <- cpg_table[,-1]

#Reading in the 450K dataset from Magnaye 2022 et al.
magnaye2022EPIC_unformatted_table <- read.table("Magnaye2022/GSE201872-GPL21145_series_matrix.txt",
                                                comment = "!",
                                                skip=5,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
magnaye2022EPIC_unformatted_samples <- magnaye2022EPIC_unformatted_table[2,-1]
magnaye2022EPIC_formatted_samples <- t(magnaye2022EPIC_unformatted_samples)[,1]

#Formatting ages.
magnaye2022EPIC_unformatted_ages <- magnaye2022EPIC_unformatted_table[3,-1]
magnaye2022EPIC_formatted_ages <- as.numeric(str_sub(magnaye2022EPIC_unformatted_ages, 6, 7))

#Formatting sex.
magnaye2022EPIC_unformatted_sex <- magnaye2022EPIC_unformatted_table[4,-1]
magnaye2022EPIC_formatted_sex <- str_sub(magnaye2022EPIC_unformatted_sex, 6, 15)

#Formatting condition.
magnaye2022EPIC_unformatted_condition <- magnaye2022EPIC_unformatted_table[5,-1]
magnaye2022EPIC_formatted_condition <- str_sub(magnaye2022EPIC_unformatted_condition, 9, 18)

#Looks like the values provided are heavily pre-processed m-values. To ensure consistency,
# I will process raw data to generate normalized beta values.
magnaye2022EPIC_list_of_files<-list.files(file.path("Magnaye2022/raw_data"))
magnaye2022EPIC_list_of_files <- magnaye2022EPIC_list_of_files[231:length(magnaye2022EPIC_list_of_files)]
magnaye2022EPIC_parsed_list_of_files <- substr(magnaye2022EPIC_list_of_files,1,30)
magnaye2022EPIC_parsed_list_of_files <- unique(magnaye2022EPIC_parsed_list_of_files)
magnaye2022EPIC_samplesheet <- data.frame(Sample=magnaye2022EPIC_formatted_samples,each=2,
                                          Ages = magnaye2022EPIC_formatted_ages,each=2,
                                          Sex= magnaye2022EPIC_formatted_sex,each=2,
                                          Condition = magnaye2022EPIC_formatted_condition,each=2,
                                          Basename = magnaye2022EPIC_parsed_list_of_files)
setwd("Magnaye2022/raw_data")
magnaye2022EPIC_RGSet <- read.metharray.exp(targets = magnaye2022EPIC_samplesheet)
magnaye2022EPIC_MSet <- preprocessSWAN(magnaye2022EPIC_RGSet)
magnaye2022EPIC_cpgs <- getBeta(magnaye2022EPIC_MSet)
setwd("..")
setwd("..")
colnames(magnaye2022EPIC_cpgs) <- magnaye2022EPIC_formatted_samples
magnaye2022EPIC_cpgs <- data.frame(magnaye2022EPIC_cpgs)

#Reading in the EPIC dataset from Magnaye 2022 et al.
magnaye2022450k_unformatted_table <- read.table("Magnaye2022/GSE201872-GPL13534_series_matrix.txt",
                                                comment = "!",
                                                skip=5,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
magnaye2022450k_unformatted_samples <- magnaye2022450k_unformatted_table[1,-1]
magnaye2022450k_formatted_samples <- t(magnaye2022450k_unformatted_samples)[,1]

#Formatting ages.
magnaye2022450k_unformatted_ages <- magnaye2022450k_unformatted_table[2,-1]
magnaye2022450k_formatted_ages <- as.numeric(str_sub(magnaye2022450k_unformatted_ages, 6, 7))

#Formatting sex.
magnaye2022450k_unformatted_sex <- magnaye2022450k_unformatted_table[3,-1]
magnaye2022450k_formatted_sex <- str_sub(magnaye2022450k_unformatted_sex, 6, 15)

#Formatting condition.
magnaye2022450k_unformatted_condition <- magnaye2022450k_unformatted_table[4,-1]
magnaye2022450k_formatted_condition <- str_sub(magnaye2022450k_unformatted_condition, 9, 18)

#Looks like the values provided are heavily pre-processed m-values. To ensure consistency,
# I will process raw data to generate normalized beta values.
magnaye2022450k_list_of_files<-list.files(file.path("Magnaye2022/raw_data"))
magnaye2022450k_list_of_files <- magnaye2022450k_list_of_files[7:230]
magnaye2022450k_parsed_list_of_files <- substr(magnaye2022450k_list_of_files,1,28)
magnaye2022450k_parsed_list_of_files <- unique(magnaye2022450k_parsed_list_of_files)
magnaye2022450k_samplesheet <- data.frame(Sample=magnaye2022450k_formatted_samples,each=2,
                                          Ages = magnaye2022450k_formatted_ages,each=2,
                                          Sex= magnaye2022450k_formatted_sex,each=2,
                                          Condition = magnaye2022450k_formatted_condition,each=2,
                                          Basename = magnaye2022450k_parsed_list_of_files)
setwd("Magnaye2022/raw_data")
magnaye2022450k_RGSet <- read.metharray.exp(targets = magnaye2022450k_samplesheet)
magnaye2022450k_MSet <- preprocessSWAN(magnaye2022450k_RGSet)
magnaye2022450k_cpgs <- getBeta(magnaye2022450k_MSet)
setwd("..")
setwd("..")
colnames(magnaye2022450k_cpgs) <- magnaye2022450k_formatted_samples
magnaye2022450k_cpgs <- data.frame(magnaye2022450k_cpgs)

magnaye2022EPIC_cpgs <- magnaye2022EPIC_cpgs[rownames(magnaye2022EPIC_cpgs) %in%
                                               rownames(magnaye2022450k_cpgs),]
magnaye2022450k_cpgs <- magnaye2022450k_cpgs[rownames(magnaye2022450k_cpgs) %in%
                                               rownames(magnaye2022EPIC_cpgs),]
shared_cpgs <- rownames(magnaye2022450k_cpgs)
#Initialize sample table.
sample_table <- data.frame(ID=character(),
                           Author=character(),
                           Year=integer(),
                           Tissue=character(),
                           CellType=character(),
                           Age=integer(),
                           Condition=character(),
                           Sex=character(),
                           DonorID=character(),
                           Misc=character())
magnayefirst_IDs <- paste(rep("F",30),1:30,sep="")
magnayefirst_samples <- data.frame(ID=magnaye2022EPIC_formatted_samples,
                                   Author=rep("Magnaye",30),
                                   Year=rep(2022,30),
                                   Tissue=rep("Bronchi",30),
                                   CellType=rep("Epithelial",30),
                                   Age=magnaye2022EPIC_formatted_ages,
                                   Condition=magnaye2022EPIC_formatted_condition,
                                   Sex=magnaye2022EPIC_formatted_sex,
                                   DonorID=magnayefirst_IDs,
                                   Misc=rep("",30))
magnayesecond_IDs <- paste(rep("G",112),1:112,sep="")
magnayesecond_samples <- data.frame(ID=magnaye2022450k_formatted_samples,
                                    Author=rep("Magnaye",112),
                                    Year=rep(2022,112),
                                    Tissue=rep("Bronchi",112),
                                    CellType=rep("Epithelial",112),
                                    Age=magnaye2022450k_formatted_ages,
                                    Condition=magnaye2022450k_formatted_condition,
                                    Sex=magnaye2022450k_formatted_sex,
                                    DonorID=magnayesecond_IDs,
                                    Misc=rep("",112))
sample_table <- rbind(magnayefirst_samples,magnayesecond_samples)
magnaye2022450k_cpgs$cpg <- rownames(magnaye2022450k_cpgs)
magnaye2022EPIC_cpgs$cpg <- rownames(magnaye2022EPIC_cpgs)
magnaye2022EPIC_cpgs <- data.table(magnaye2022EPIC_cpgs)
magnaye2022450k_cpgs <- data.table(magnaye2022450k_cpgs)
cpg_table <- merge(magnaye2022450k_cpgs,magnaye2022EPIC_cpgs)

# Time for a monocyte data set.
estupinan2022_unformatted_table <- read.table("Estupinan2022/GSE201752_series_matrix.txt",
                                              comment = "!",
                                              skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
estupinan2022_formatted_samples <- strsplit(estupinan2022_unformatted_table[1,2]," ")[[1]]

#Formatting ages.
estupinan2022_unformatted_ages <- estupinan2022_unformatted_table[6,-1]
estupinan2022_formatted_ages <- as.numeric(str_sub(estupinan2022_unformatted_ages, 6, 7))

#Formatting sex.
estupinan2022_unformatted_sex <- estupinan2022_unformatted_table[5,-1]
estupinan2022_formatted_sex <- str_sub(estupinan2022_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
estupinan2022_unformatted_condition <- estupinan2022_unformatted_table[4,-1]
estupinan2022_formatted_condition <- str_sub(estupinan2022_unformatted_condition, 17, 20)
estupinan2022_formatted_condition[estupinan2022_formatted_condition == "HD"] <- "Control"
estupinan2022_formatted_misc <- estupinan2022_formatted_condition
estupinan2022_formatted_condition[estupinan2022_formatted_condition != "Control"] <- "Giant Cell Arteritis"
estupinan2022_IDs <- paste(rep("H",113),1:113,sep="")

estupinan2022_samples <- data.frame(ID=estupinan2022_formatted_samples,
                                    Author=rep("Estupinan",113),
                                    Year=rep(2022,113),
                                    Tissue=rep("Blood",113),
                                    CellType=rep("Monocytes",113),
                                    Age=estupinan2022_formatted_ages,
                                    Condition=estupinan2022_formatted_condition,
                                    Sex=estupinan2022_formatted_sex,
                                    DonorID=estupinan2022_IDs,
                                    Misc=estupinan2022_formatted_misc)

estupinan2022_cpgs <- data.frame(read_table2("Estupinan2022/GSE201752_processed_data.txt", skip=4),row.names=1)
colnames(estupinan2022_cpgs) <- estupinan2022_formatted_samples
estupinan2022_cpgs <- estupinan2022_cpgs[rownames(estupinan2022_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(estupinan2022_cpgs),]

estupinan2022_cpgs$cpg <- rownames(estupinan2022_cpgs)
estupinan2022_cpgs <- data.table(estupinan2022_cpgs)
sample_table <- rbind(sample_table,estupinan2022_samples)
cpg_table <- merge(cpg_table,estupinan2022_cpgs, all=TRUE)
#
#
# #Time for a blood data set.
okereke2021_unformatted_table <- read.table("Okereke2021/GSE190540_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
okereke2021_unformatted_samples <- strsplit(okereke2021_unformatted_table[1,2]," ")[[1]]

#Formatting ages.
okereke2021_unformatted_ages <- okereke2021_unformatted_table[4,-1]
okereke2021_formatted_ages <- as.numeric(str_sub(okereke2021_unformatted_ages, 6, 7))

#Formatting sex.
okereke2021_unformatted_sex <- okereke2021_unformatted_table[2,-1]
okereke2021_formatted_sex <- str_sub(okereke2021_unformatted_sex, 9, 15)

#Formatting condition and miscellaneous information.
okereke2021_unformatted_condition <- okereke2021_unformatted_table[3,-1]
okereke2021_formatted_condition <- str_sub(okereke2021_unformatted_condition, 17, 25)
okereke2021_formatted_condition[okereke2021_formatted_condition == "Case"] <- "Cognitive Impairment"

okereke2021_IDs <- paste(rep("I",90),1:90,sep="")

okereke2021_samples <- data.frame(ID=okereke2021_unformatted_samples,
                                  Author=rep("Okereke",90),
                                  Year=rep(2021,90),
                                  Tissue=rep("Blood",90),
                                  CellType=rep("PBMCs",90),
                                  Age=okereke2021_formatted_ages,
                                  Condition=okereke2021_formatted_condition,
                                  Sex=okereke2021_formatted_sex,
                                  DonorID=okereke2021_IDs,
                                  Misc=rep(NA,90))

okereke2021_cpgs <- read.table("Okereke2021/GSE190540_series_matrix.txt",
                               comment = "!",
                               skip=5,
                               fill=TRUE)
okereke2021_cpgs <- data.frame(okereke2021_cpgs)
okereke2021_cpgs <- okereke2021_cpgs[-c(1:5),]
rownames(okereke2021_cpgs) <- okereke2021_cpgs$V1
okereke2021_cpgs <- data.frame(okereke2021_cpgs[,-1])
colnames(okereke2021_cpgs) <- okereke2021_unformatted_samples
okereke2021_cpgs <- okereke2021_cpgs[rownames(okereke2021_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(okereke2021_cpgs),]

okereke2021_cpgs$cpg <- rownames(okereke2021_cpgs)
okereke2021_cpgs <- data.table(okereke2021_cpgs)
sample_table <- rbind(sample_table,okereke2021_samples)
cpg_table <- merge(cpg_table,okereke2021_cpgs, all=TRUE)
#
# # Skin data set.
muse2021_unformatted_table <- read.table("Muse2021/GSE188593_series_matrix.txt",
                                         comment = "!",
                                         skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
muse2021_unformatted_samples <- muse2021_unformatted_table[1,-1]
muse2021_formatted_samples <- t(muse2021_unformatted_samples)[,1]

#Formatting ages.
muse2021_unformatted_ages <- muse2021_unformatted_table[3,-1]
muse2021_formatted_ages <- as.numeric(str_sub(muse2021_unformatted_ages, 13, 19))

#Formatting sex.
muse2021_unformatted_sex <- muse2021_unformatted_table[4,-1]
muse2021_formatted_sex <- str_sub(muse2021_unformatted_sex, 6, 7)

muse2021_unformatted_donor <- muse2021_unformatted_table[2,-1]
muse2021_unformatted_donor <- str_sub(muse2021_unformatted_donor, 10,15)
muse2021_unformatted_donor <- as.numeric(factor(muse2021_unformatted_donor))
muse2021_formatted_donor <- paste(rep("J",64),muse2021_unformatted_donor,sep="")

#Formatting condition and miscellaneous information.
muse2021_unformatted_condition <- muse2021_unformatted_table[6,-1]
muse2021_formatted_condition <- str_sub(muse2021_unformatted_condition, 17, 30)
muse2021_formatted_condition[muse2021_formatted_condition == ""] <- "Control"

muse2021_IDs <- paste(rep("J",64),1:64,sep="")

muse2021_samples <- data.frame(ID=muse2021_formatted_samples,
                               Author=rep("Muse",64),
                               Year=rep(2021,64),
                               Tissue=rep("Skin",64),
                               CellType=rep("Epithelial",64),
                               Age=muse2021_formatted_ages,
                               Condition=muse2021_formatted_condition,
                               Sex=muse2021_formatted_sex,
                               DonorID=muse2021_IDs,
                               Misc=rep(NA,64))

muse2021_cpgs <- read.table("Muse2021/GSE188593_series_matrix.txt",
                            comment = "!",
                            skip=5,
                            fill=TRUE)
muse2021_cpgs <- data.frame(muse2021_cpgs)
muse2021_cpgs <- muse2021_cpgs[-c(1:6),]
rownames(muse2021_cpgs) <- muse2021_cpgs$V1
muse2021_cpgs <- data.frame(muse2021_cpgs[,-1])
colnames(muse2021_cpgs) <- muse2021_unformatted_samples
muse2021_cpgs <- muse2021_cpgs[rownames(muse2021_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(muse2021_cpgs),]

muse2021_cpgs$cpg <- rownames(muse2021_cpgs)
muse2021_cpgs <- data.table(muse2021_cpgs)
sample_table <- rbind(sample_table,muse2021_samples)
cpg_table <- merge(cpg_table,muse2021_cpgs, all=TRUE)

# #Blood data set.
konigsberg2021_unformatted_table <- read.table("Konigsberg2021/GSE167202_series_matrix.txt",
                                               comment = "!",
                                               skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
konigsberg2021_formatted_samples <- strsplit(konigsberg2021_unformatted_table[1,2]," ")[[1]]

#Formatting ages.
konigsberg2021_unformatted_ages <- konigsberg2021_unformatted_table[4,-1]
konigsberg2021_formatted_ages <- as.numeric(str_sub(konigsberg2021_unformatted_ages, 5, 8))

#Formatting sex.
konigsberg2021_unformatted_sex <- konigsberg2021_unformatted_table[3,-1]
konigsberg2021_formatted_sex <- str_sub(konigsberg2021_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
konigsberg2021_unformatted_condition <- konigsberg2021_unformatted_table[2,-1]
konigsberg2021_formatted_condition <- str_sub(konigsberg2021_unformatted_condition, 15, 30)
konigsberg2021_formatted_condition[konigsberg2021_formatted_condition == "negative"] <- "Control"
konigsberg2021_formatted_condition[konigsberg2021_formatted_condition == "positive"] <- "COVID"
konigsberg2021_formatted_condition[konigsberg2021_formatted_condition == "other infection"] <-
  "Respiratory Illness"

konigsberg2021_IDs <- paste(rep("K",525),1:525,sep="")

konigsberg2021_samples <- data.frame(ID=konigsberg2021_formatted_samples,
                                     Author=rep("Konigsberg",525),
                                     Year=rep(2021,525),
                                     Tissue=rep("Blood",525),
                                     CellType=rep("PBMCs",525),
                                     Age=konigsberg2021_formatted_ages,
                                     Condition=konigsberg2021_formatted_condition,
                                     Sex=konigsberg2021_formatted_sex,
                                     DonorID=konigsberg2021_IDs,
                                     Misc=rep(NA,525))

konigsberg2021_cpgs <- data.table::fread("Konigsberg2021/GSE167202_ProcessedBetaValues.txt",
                                         header=TRUE) %>%
  as.data.frame()
konigsberg2021_cpgs <- konigsberg2021_cpgs[,-c(2:22)]
konigsberg2021_sequence_ids <- konigsberg2021_unformatted_table[5,-1]
konigsberg2021_sequence_ids <- unlist(konigsberg2021_sequence_ids)
row.names(konigsberg2021_cpgs) <- konigsberg2021_cpgs$ID_REF
konigsberg2021_cpgs <- konigsberg2021_cpgs[,-1]
konigsberg2021_cpgs <- konigsberg2021_cpgs[,konigsberg2021_sequence_ids]
colnames(konigsberg2021_cpgs) <- konigsberg2021_formatted_samples

konigsberg2021_cpgs <- konigsberg2021_cpgs[rownames(konigsberg2021_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(konigsberg2021_cpgs),]

konigsberg2021_cpgs$cpg <- rownames(konigsberg2021_cpgs)
konigsberg2021_cpgs <- data.table(konigsberg2021_cpgs)
sample_table <- rbind(sample_table,konigsberg2021_samples)
cpg_table <- merge(cpg_table,konigsberg2021_cpgs, all=TRUE)

#Brain samples.
haghighi2022_unformatted_table <- read.table("Haghighi2022/GSE191200_series_matrix.txt",
                                             comment = "!",
                                             skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
haghighi2022_formatted_samples <- strsplit(haghighi2022_unformatted_table[1,2]," ")[[1]]

#Formatting ages.
haghighi2022_unformatted_ages <- haghighi2022_unformatted_table[4,-1]
haghighi2022_formatted_ages <- as.numeric(str_sub(haghighi2022_unformatted_ages, 5, 8))

#Formatting sex.
haghighi2022_unformatted_sex <- haghighi2022_unformatted_table[5,-1]
haghighi2022_formatted_sex <- str_sub(haghighi2022_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
haghighi2022_unformatted_condition <- haghighi2022_unformatted_table[6,-1]
haghighi2022_formatted_condition <- str_sub(haghighi2022_unformatted_condition, 12, 30)

haghighi2022_unformatted_donor <- haghighi2022_unformatted_table[2,-1]
haghighi2022_unformatted_donor <- str_sub(haghighi2022_unformatted_donor,11,15)
haghighi2022_unformatted_donor <- as.numeric(factor(haghighi2022_unformatted_donor,
                                                    labels=c(1:22)))
haghighi2022_formatted_donor <- paste(rep("L",56),haghighi2022_unformatted_donor,sep="")

haghighi2022_samples <- data.frame(ID=haghighi2022_formatted_samples,
                                   Author=rep("Haghighi",56),
                                   Year=rep(2022,56),
                                   Tissue=rep("Brain",56),
                                   CellType=rep("Microglia",56),
                                   Age=haghighi2022_formatted_ages,
                                   Condition=haghighi2022_formatted_condition,
                                   Sex=haghighi2022_formatted_sex,
                                   DonorID=haghighi2022_formatted_donor,
                                   Misc=rep(NA,56))

haghighi2022_cpgs <- read.table("Haghighi2022/GSE191200_series_matrix.txt",
                                comment = "!",
                                skip=7,
                                fill=TRUE)
haghighi2022_cpgs <- haghighi2022_cpgs[-c(1:6),]
row.names(haghighi2022_cpgs) <- haghighi2022_cpgs$V1
haghighi2022_cpgs <- haghighi2022_cpgs[,-1]
colnames(haghighi2022_cpgs) <- haghighi2022_cpgs[1,]
haghighi2022_cpgs <- haghighi2022_cpgs[-1,]

haghighi2022_cpgs <- haghighi2022_cpgs[rownames(haghighi2022_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(haghighi2022_cpgs),]

haghighi2022_cpgs$cpg <- rownames(haghighi2022_cpgs)
haghighi2022_cpgs <- data.table(haghighi2022_cpgs)
sample_table <- rbind(sample_table,haghighi2022_samples)
cpg_table <- merge(cpg_table,haghighi2022_cpgs, all=TRUE)

# #Colorectal samples.
chen2021_unformatted_table <- read.table("Chen2021/GSE159898_series_matrix.txt",
                                         comment = "!",
                                         skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
chen2021_formatted_samples <- strsplit(chen2021_unformatted_table[1,2]," ")[[1]]

#Formatting ages.
chen2021_unformatted_ages <- chen2021_unformatted_table[3,-1]
chen2021_formatted_ages <- as.numeric(str_sub(chen2021_unformatted_ages, 5, 8))

#Formatting sex.
chen2021_unformatted_sex <- chen2021_unformatted_table[2,-1]
chen2021_formatted_sex <- str_sub(chen2021_unformatted_sex, 9, 15)

#Formatting condition and miscellaneous information.
chen2021_unformatted_condition <- chen2021_unformatted_table[4,-1]
chen2021_formatted_condition <- str_sub(chen2021_unformatted_condition, 14, 40)
chen2021_formatted_condition[chen2021_formatted_condition == "normal colorectal tissue"] <-
  "Control"
chen2021_formatted_condition[chen2021_formatted_condition == "colorectal cancer tissue"] <-
  "Colorectal Cancer"

chen2021_unformatted_donor <- chen2021_unformatted_table[5,-1]
chen2021_unformatted_donor <- str_sub(chen2021_unformatted_donor,11,15)
chen2021_unformatted_donor <- as.numeric(factor(chen2021_unformatted_donor,
                                                labels=c(1:21)))
chen2021_formatted_donor <- paste(rep("M",44),chen2021_unformatted_donor,sep="")

chen2021_samples <- data.frame(ID=chen2021_formatted_samples,
                               Author=rep("Chen",44),
                               Year=rep(2021,44),
                               Tissue=rep("Colorectal",44),
                               CellType=rep("Epithelial (Colorectal)",44),
                               Age=chen2021_formatted_ages,
                               Condition=chen2021_formatted_condition,
                               Sex=chen2021_formatted_sex,
                               DonorID=chen2021_formatted_donor,
                               Misc=rep(NA,44))

chen2021_cpgs <- read.table("Chen2021/GSE159898_series_matrix.txt",
                            comment = "!",
                            skip=7,
                            fill=TRUE)
chen2021_cpgs <- chen2021_cpgs[-c(1:5),]
row.names(chen2021_cpgs) <- chen2021_cpgs$V1
chen2021_cpgs <- chen2021_cpgs[,-1]
colnames(chen2021_cpgs) <- chen2021_cpgs[1,]
chen2021_cpgs <- chen2021_cpgs[-1,]

chen2021_cpgs <- chen2021_cpgs[rownames(chen2021_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(chen2021_cpgs),]

chen2021_cpgs$cpg <- rownames(chen2021_cpgs)
chen2021_cpgs <- data.table(chen2021_cpgs)
sample_table <- rbind(sample_table,chen2021_samples)
cpg_table <- merge(cpg_table,chen2021_cpgs, all=TRUE)

#Skeletal Muscle samples
voisin2021_unformatted_table <- read.table("Voisin2021/GSE151407_series_matrix.txt",
                                           comment = "!",
                                           skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
voisin2021_formatted_samples <- strsplit(voisin2021_unformatted_table[1,2]," ")[[1]]

#Formatting ages.
voisin2021_unformatted_ages <- voisin2021_unformatted_table[3,-1]
voisin2021_formatted_ages <- as.numeric(str_sub(voisin2021_unformatted_ages, 5, 8))

#Formatting sex.
voisin2021_unformatted_sex <- voisin2021_unformatted_table[4,-1]
voisin2021_formatted_sex <- str_sub(voisin2021_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
voisin2021_unformatted_condition <- voisin2021_unformatted_table[2,-1]
voisin2021_formatted_condition <- str_sub(voisin2021_unformatted_condition, 8, 10)
voisin2021_formatted_condition[voisin2021_formatted_condition == "PRE"] <-
  "Control"
voisin2021_formatted_condition[voisin2021_formatted_condition != "Control"] <-
  "HIIT"

voisin2021_unformatted_donor <- voisin2021_unformatted_table[2,-1]
voisin2021_unformatted_donor <- str_sub(voisin2021_unformatted_donor,-9,-6)
voisin2021_unformatted_donor <- sub("_", "", voisin2021_unformatted_donor)
voisin2021_unformatted_donor <- as.numeric(factor(voisin2021_unformatted_donor,
                                                  labels=c(1:25)))
voisin2021_formatted_donor <- paste(rep("N",78),voisin2021_unformatted_donor,sep="")

voisin2021_samples <- data.frame(ID=voisin2021_formatted_samples,
                                 Author=rep("Voisin",78),
                                 Year=rep(2021,78),
                                 Tissue=rep("Skeletal Muscle",78),
                                 CellType=rep("Muscle Cells",78),
                                 Age=voisin2021_formatted_ages,
                                 Condition=voisin2021_formatted_condition,
                                 Sex=voisin2021_formatted_sex,
                                 DonorID=voisin2021_formatted_donor,
                                 Misc=rep(NA,78))

voisin2021_cpgs <- read.table("Voisin2021/GSE151407_series_matrix.txt",
                              comment = "!",
                              skip=7,
                              fill=TRUE)
voisin2021_cpgs <- voisin2021_cpgs[-c(1:4),]
row.names(voisin2021_cpgs) <- voisin2021_cpgs$V1
voisin2021_cpgs <- voisin2021_cpgs[,-1]
colnames(voisin2021_cpgs) <- voisin2021_cpgs[1,]
voisin2021_cpgs <- voisin2021_cpgs[-1,]

voisin2021_cpgs <- voisin2021_cpgs[rownames(voisin2021_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(voisin2021_cpgs),]

voisin2021_cpgs$cpg <- rownames(voisin2021_cpgs)
voisin2021_cpgs <- data.table(voisin2021_cpgs)
sample_table <- rbind(sample_table,voisin2021_samples)
cpg_table <- merge(cpg_table,voisin2021_cpgs, all=TRUE)

# More brain. Note - I am somewhat apprehensive about this following dataset.

fries2019_unformatted_table <- read.table("Fries2019/GSE129428_series_matrix.txt",
                                          comment = "!",
                                          skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
fries2019_formatted_samples <- strsplit(fries2019_unformatted_table[1,2]," ")[[1]]

#Formatting ages.
fries2019_unformatted_ages <- fries2019_unformatted_table[4,-1]
fries2019_formatted_ages <- as.numeric(str_sub(fries2019_unformatted_ages, 14, 19))

#Formatting sex.
fries2019_unformatted_sex <- fries2019_unformatted_table[2,-1]
fries2019_formatted_sex <- str_sub(fries2019_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
fries2019_unformatted_condition <- fries2019_unformatted_table[3,-1]
fries2019_formatted_condition <- str_sub(fries2019_unformatted_condition, 8, 100)

fries2019_formatted_donor <- paste(rep("O",64),c(1:64),sep="")

fries2019_samples <- data.frame(ID=fries2019_formatted_samples,
                                Author=rep("Fries",64),
                                Year=rep(2019,64),
                                Tissue=rep("Brain",64),
                                CellType=rep("Brain Cells",64),
                                Age=fries2019_formatted_ages,
                                Condition=fries2019_formatted_condition,
                                Sex=fries2019_formatted_sex,
                                DonorID=fries2019_formatted_donor,
                                Misc=rep(NA,64))

fries2019_cpgs <- read.table("Fries2019/GSE129428_series_matrix.txt",
                             comment = "!",
                             skip=7,
                             fill=TRUE)
fries2019_cpgs <- fries2019_cpgs[-c(1:4),]
row.names(fries2019_cpgs) <- fries2019_cpgs$V1
fries2019_cpgs <- fries2019_cpgs[,-1]
colnames(fries2019_cpgs) <- fries2019_cpgs[1,]
fries2019_cpgs <- fries2019_cpgs[-1,]

fries2019_cpgs <- fries2019_cpgs[rownames(fries2019_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(fries2019_cpgs),]

fries2019_cpgs$cpg <- rownames(fries2019_cpgs)
fries2019_cpgs <- data.table(fries2019_cpgs)
sample_table <- rbind(sample_table,fries2019_samples)
cpg_table <- merge(cpg_table,fries2019_cpgs, all=TRUE)

#Samples from children.

islam2018_unformatted_table <- read.table("Islam2018/GSE124366_series_matrix.txt",
                                          comment = "!",
                                          skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
islam2018_unformatted_samples <- islam2018_unformatted_table[6,]
islam2018_formatted_samples <- str_sub(islam2018_unformatted_table[6,],1,100)[2:216]

#Formatting ages.
islam2018_unformatted_ages <- islam2018_unformatted_table[4,-1]
islam2018_formatted_ages <- as.numeric(str_sub(islam2018_unformatted_ages, 35,100))

#Formatting sex.
islam2018_unformatted_sex <- islam2018_unformatted_table[3,-1]
islam2018_formatted_sex <- str_sub(islam2018_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
islam2018_formatted_condition <- rep("Control",215)

islam2018_unformatted_donor <- islam2018_unformatted_table[1,-1]
islam2018_unformatted_donor <- str_sub(islam2018_unformatted_donor,1,8)
islam2018_unformatted_donor <- as.numeric(factor(islam2018_unformatted_donor,
                                                 labels=c(1:201)))
islam2018_formatted_donor <- paste(rep("P",201),islam2018_unformatted_donor,sep="")

islam2018_unformatted_celltype <- islam2018_unformatted_table[2,-1]
islam2018_formatted_celltype <- str_sub(islam2018_unformatted_celltype, 9, 20)

islam2018_unformatted_tissue <- islam2018_formatted_celltype
islam2018_unformatted_tissue[islam2018_unformatted_tissue == "PBMC"] <- "Blood"
islam2018_formatted_tissue <- islam2018_unformatted_tissue

islam2018_samples <- data.frame(ID=islam2018_formatted_samples,
                                Author=rep("Islam",215),
                                Year=rep(2018,215),
                                Tissue=islam2018_formatted_tissue,
                                CellType=islam2018_formatted_celltype,
                                Age=islam2018_formatted_ages,
                                Condition=islam2018_formatted_condition,
                                Sex=islam2018_formatted_sex,
                                DonorID=islam2018_formatted_donor,
                                Misc=rep(NA,215))

islam2018_cpgs <- read.table("Islam2018/GSE124366_series_matrix.txt",
                             comment = "!",
                             skip=6,
                             fill=TRUE)
islam2018_cpgs <- islam2018_cpgs[-c(1:6),]
row.names(islam2018_cpgs) <- islam2018_cpgs$V1
islam2018_cpgs <- islam2018_cpgs[,-1]
colnames(islam2018_cpgs) <- islam2018_formatted_samples

islam2018_cpgs <- islam2018_cpgs[rownames(islam2018_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(islam2018_cpgs),]

islam2018_cpgs$cpg <- rownames(islam2018_cpgs)
islam2018_cpgs <- data.table(islam2018_cpgs)
sample_table <- rbind(sample_table,islam2018_samples)
cpg_table <- merge(cpg_table,islam2018_cpgs, all=TRUE)
#

#Samples from breast tissue.

johnson2016_unformatted_table <- read.table("Johnson2016/GSE88883_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 4)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
johnson2016_unformatted_samples <- johnson2016_unformatted_table[4,]
johnson2016_formatted_samples <- str_sub(johnson2016_unformatted_samples,1,100)[2:101]

#Formatting ages.
johnson2016_unformatted_ages <- johnson2016_unformatted_table[2,-1]
johnson2016_formatted_ages <- as.numeric(str_sub(johnson2016_unformatted_ages, 14,100))

#Formatting sex.
johnson2016_unformatted_sex <- johnson2016_unformatted_table[1,-1]
johnson2016_formatted_sex <- str_sub(johnson2016_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
johnson2016_formatted_condition <- rep("Control",100)

johnson2016_formatted_donor <- paste(rep("Q",100),c(1:100),sep="")

johnson2016_samples <- data.frame(ID=johnson2016_formatted_samples,
                                  Author=rep("Johnson",100),
                                  Year=rep(2016,100),
                                  Tissue=rep("Breast",100),
                                  CellType=rep("Breast",100),
                                  Age=johnson2016_formatted_ages,
                                  Condition=johnson2016_formatted_condition,
                                  Sex=johnson2016_formatted_sex,
                                  DonorID=johnson2016_formatted_donor,
                                  Misc=rep(NA,100))

johnson2016_cpgs <- read.table("Johnson2016/GSE88883_series_matrix.txt",
                               comment = "!",
                               skip=4,
                               fill=TRUE)
johnson2016_cpgs <- johnson2016_cpgs[-c(1:4),]
row.names(johnson2016_cpgs) <- johnson2016_cpgs$V1
johnson2016_cpgs <- johnson2016_cpgs[,-1]
colnames(johnson2016_cpgs) <- johnson2016_formatted_samples

johnson2016_cpgs <- johnson2016_cpgs[rownames(johnson2016_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(johnson2016_cpgs),]

johnson2016_cpgs$cpg <- rownames(johnson2016_cpgs)
johnson2016_cpgs <- data.table(johnson2016_cpgs)
sample_table <- rbind(sample_table,johnson2016_samples)
cpg_table <- merge(cpg_table,johnson2016_cpgs, all=TRUE)

#Samples from liver.

horvath2014_unformatted_table <- read.table("Horvath2014/GSE61258_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 4)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
horvath2014_unformatted_samples <- horvath2014_unformatted_table[4,]
horvath2014_formatted_samples <- str_sub(horvath2014_unformatted_samples,1,100)[2:80]

#Formatting ages.
horvath2014_unformatted_ages <- horvath2014_unformatted_table[2,-1]
horvath2014_formatted_ages <- as.numeric(str_sub(horvath2014_unformatted_ages, 5,100))

#Formatting sex.
horvath2014_unformatted_sex <- horvath2014_unformatted_table[1,-1]
horvath2014_formatted_sex <- str_sub(horvath2014_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
horvath2014_formatted_condition <- rep("Control",79)

horvath2014_formatted_donor <- paste(rep("R",79),c(1:79),sep="")

horvath2014_samples <- data.frame(ID=horvath2014_formatted_samples,
                                  Author=rep("Horvath",79),
                                  Year=rep(2014,79),
                                  Tissue=rep("Liver",79),
                                  CellType=rep("Liver",79),
                                  Age=horvath2014_formatted_ages,
                                  Condition=horvath2014_formatted_condition,
                                  Sex=horvath2014_formatted_sex,
                                  DonorID=horvath2014_formatted_donor,
                                  Misc=rep(NA,79))

horvath2014_cpgs <- read.table("Horvath2014/GSE61258_series_matrix.txt",
                               comment = "!",
                               skip=4,
                               fill=TRUE)
horvath2014_cpgs <- horvath2014_cpgs[-c(1:4),]
row.names(horvath2014_cpgs) <- horvath2014_cpgs$V1
horvath2014_cpgs <- horvath2014_cpgs[,-1]
colnames(horvath2014_cpgs) <- horvath2014_formatted_samples

horvath2014_cpgs <- horvath2014_cpgs[rownames(horvath2014_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(horvath2014_cpgs),]

horvath2014_cpgs$cpg <- rownames(horvath2014_cpgs)
horvath2014_cpgs <- data.table(horvath2014_cpgs)
sample_table <- rbind(sample_table,horvath2014_samples)
cpg_table <- merge(cpg_table,horvath2014_cpgs, all=TRUE)

#Samples from sperm tissue.

pilsner2022_unformatted_table <- read.table("Pilsner2022/GSE185445_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
pilsner2022_unformatted_samples <- pilsner2022_unformatted_table[1,-1]
pilsner2022_formatted_samples <- str_sub(pilsner2022_unformatted_samples,1,100)[1:379]

#Formatting ages.
pilsner2022_unformatted_ages <- pilsner2022_unformatted_table[3,-1]
pilsner2022_formatted_ages <- as.numeric(str_sub(pilsner2022_unformatted_ages, 8,100))

#Formatting sex.
pilsner2022_formatted_sex <- rep("Male",379)

#Formatting condition and miscellaneous information.
pilsner2022_formatted_condition <- rep("Control",379)

pilsner2022_formatted_donor <- paste(rep("T",379),c(1:379),sep="")

pilsner2022_samples <- data.frame(ID=pilsner2022_formatted_samples,
                                  Author=rep("Pilsner",379),
                                  Year=rep(2022,379),
                                  Tissue=rep("Semen",379),
                                  CellType=rep("Sperm",379),
                                  Age=pilsner2022_formatted_ages,
                                  Condition=pilsner2022_formatted_condition,
                                  Sex=pilsner2022_formatted_sex,
                                  DonorID=pilsner2022_formatted_donor,
                                  Misc=rep(NA,379))

pilsner2022_unformatted_basenames <- pilsner2022_unformatted_table[2,-1]
pilsner2022_formatted_basenames <- str_sub(pilsner2022_unformatted_basenames,11,100)

pilsner2022_cpgs <- data.table::fread("Pilsner2022/GSE185445_processed_methylation_data_matrix.txt",
                                      header=FALSE)  %>%
  as.data.frame()

pilsner2022_colnames <- data.frame(read_table2("Pilsner2022/colnames.csv",
                                               col_names = FALSE))
pilsner2022_colnames <- substr(pilsner2022_colnames,2,20)

row.names(pilsner2022_cpgs) <- pilsner2022_cpgs[,1]
pilsner2022_cpgs <- pilsner2022_cpgs[,-1]
colnames(pilsner2022_cpgs) <- pilsner2022_colnames
pilsner2022_cpgs <- pilsner2022_cpgs[,pilsner2022_formatted_basenames]

pilsner2022_cpgs <- pilsner2022_cpgs[rownames(pilsner2022_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(pilsner2022_cpgs),]
colnames(pilsner2022_cpgs) <- pilsner2022_formatted_samples

pilsner2022_cpgs$cpg <- rownames(pilsner2022_cpgs)
pilsner2022_cpgs <- data.table(pilsner2022_cpgs)
sample_table <- rbind(sample_table,pilsner2022_samples)
cpg_table <- merge(cpg_table,pilsner2022_cpgs, all=TRUE)

#Samples from nasal epithelia

davalos2022_unformatted_table <- read.table("Davalos2022/GSE193879_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 4)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
davalos2022_unformatted_samples <- davalos2022_unformatted_table[1,]
davalos2022_formatted_samples <- str_sub(davalos2022_unformatted_samples,1,62)[2:128]

#Formatting ages.
davalos2022_unformatted_ages <- davalos2022_unformatted_table[4,-1]
davalos2022_formatted_ages <- as.numeric(str_sub(davalos2022_unformatted_ages, 13,15))

#Formatting sex.
davalos2022_unformatted_sex <- davalos2022_unformatted_table[3,-1]
davalos2022_formatted_sex <- str_sub(davalos2022_unformatted_sex, 17, 30)

#Formatting condition and miscellaneous information.
davalos2022_unformatted_condition <- davalos2022_unformatted_table[2,-1]
davalos2022_formatted_condition <- str_sub(davalos2022_unformatted_condition, 8,39)

davalos2022_formatted_donor <- paste(rep("U",127),c(1:127),sep="")

davalos2022_samples <- data.frame(ID=davalos2022_formatted_samples,
                                  Author=rep("Davalos",127),
                                  Year=rep(2022,127),
                                  Tissue=rep("Blood",127),
                                  CellType=rep("Blood",127),
                                  Age=davalos2022_formatted_ages,
                                  Condition=davalos2022_formatted_condition,
                                  Sex=davalos2022_formatted_sex,
                                  DonorID=davalos2022_formatted_donor,
                                  Misc=rep(NA,127))

davalos2022_samples <- davalos2022_samples[!is.na(davalos2022_samples$Age),]

davalos2022_cpgs <- data.table::fread("Davalos2022/GSE193879_Matrix_processed.csv",
                                      header=FALSE)  %>%
  as.data.frame()
row.names(davalos2022_cpgs) <- davalos2022_cpgs$V1
davalos2022_cpgs <- davalos2022_cpgs[-1,-1]
davalos2022_cpgs <- davalos2022_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),127)]
colnames(davalos2022_cpgs) <- davalos2022_formatted_samples
davalos2022_cpgs <- davalos2022_cpgs[,colnames(davalos2022_cpgs) %in% davalos2022_samples$ID]

davalos2022_cpgs <- davalos2022_cpgs[rownames(davalos2022_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(davalos2022_cpgs),]

davalos2022_cpgs$cpg <- rownames(davalos2022_cpgs)
davalos2022_cpgs <- data.table(davalos2022_cpgs)
sample_table <- rbind(sample_table,davalos2022_samples)
cpg_table <- merge(cpg_table,davalos2022_cpgs, all=TRUE)


hannon2021_unformatted_table <- read.table("Hannon2021/GSE152026_series_matrix.txt",
                                           comment = "!",
                                           skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
hannon2021_unformatted_samples <- hannon2021_unformatted_table[2,]
hannon2021_formatted_samples <- str_sub(hannon2021_unformatted_samples,1,62)[2:935]

#Formatting ages.
hannon2021_unformatted_ages <- hannon2021_unformatted_table[5,-1]
hannon2021_formatted_ages <- as.numeric(str_sub(hannon2021_unformatted_ages, 5,100))

#Formatting sex.
hannon2021_unformatted_sex <- hannon2021_unformatted_table[4,-1]
hannon2021_formatted_sex <- str_sub(hannon2021_unformatted_sex, 6, 30)

#Formatting condition and miscellaneous information.
hannon2021_unformatted_condition <- hannon2021_unformatted_table[3,-1]
hannon2021_formatted_condition <- str_sub(hannon2021_unformatted_condition, 12,39)
hannon2021_formatted_condition[hannon2021_formatted_condition=="Case"] <- "Schizophrenia"

hannon2021_formatted_donor <- paste(rep("V",934),c(1:934),sep="")

hannon2021_samples <- data.frame(ID=hannon2021_formatted_samples,
                                 Author=rep("Hannon",934),
                                 Year=rep(2021,934),
                                 Tissue=rep("Blood",934),
                                 CellType=rep("Blood",934),
                                 Age=hannon2021_formatted_ages,
                                 Condition=hannon2021_formatted_condition,
                                 Sex=hannon2021_formatted_sex,
                                 DonorID=hannon2021_formatted_donor,
                                 Misc=rep(NA,934))

hannon2021_samples <- hannon2021_samples[!is.na(hannon2021_samples$Age),]
hannon2021_sample_names <- str_sub(hannon2021_unformatted_table[1,-1],1,19)

hannon2021_cpgs <- data.table::fread("Hannon2021/GSE152026_EUGEI_processed_signals.csv",
                                     header=FALSE)  %>%
  as.data.frame()
row.names(hannon2021_cpgs) <- hannon2021_cpgs$V1
hannon2021_cpgs <- hannon2021_cpgs[-1,-1]
hannon2021_cpgs <- hannon2021_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),934)]
hannon2021_names <- read_csv("Hannon2021/target.txt",
                             col_names = FALSE)
hannon2021_names <- str_sub(unlist(hannon2021_names),1,500)
hannon2021_names <- hannon2021_names[2:1869]
hannon2021_names <- hannon2021_names[rep(c(rep(TRUE, 2- 1), FALSE),934)]
colnames(hannon2021_cpgs) <- hannon2021_names

hannon2021_cpgs <- hannon2021_cpgs[,hannon2021_sample_names]
colnames(hannon2021_cpgs) <- hannon2021_formatted_samples
hannon2021_cpgs <- hannon2021_cpgs[,colnames(hannon2021_cpgs) %in%
                                     hannon2021_samples$ID]
hannon2021_cpgs <- hannon2021_cpgs[rownames(hannon2021_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(hannon2021_changed_cpgs),]
hannon2021_cpgs$cpg <- rownames(hannon2021_cpgs)
hannon2021_cpgs <- data.table(hannon2021_cpgs)
sample_table <- rbind(sample_table,hannon2021_samples)
cpg_table <- merge(cpg_table,hannon2021_cpgs, all=TRUE)

martino2018_unformatted_table <- read.table("Martino2018/GSE114135-GPL23976_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
martino2018_unformatted_samples <- martino2018_unformatted_table[5,-1]
martino2018_formatted_samples <- str_sub(martino2018_unformatted_samples,1,62)[1:205]

#Formatting ages.
martino2018_unformatted_ages <- martino2018_unformatted_table[2,-1]
martino2018_formatted_ages <- as.numeric(str_sub(martino2018_unformatted_ages, 5,100))

#Formatting sex.
martino2018_unformatted_sex <- martino2018_unformatted_table[3,-1]
martino2018_formatted_sex <- str_sub(martino2018_unformatted_sex, 6, 30)

#Formatting condition and miscellaneous information.
martino2018_unformatted_condition <- martino2018_unformatted_table[4,-1]
martino2018_formatted_condition <- str_sub(martino2018_unformatted_condition, -1,-1)
martino2018_formatted_condition[martino2018_formatted_condition=="1"] <- "Stimulated"
martino2018_formatted_condition[martino2018_formatted_condition=="0"] <- "Control"

martino2018_unformatted_donor <- martino2018_unformatted_table[1,-1]
martino2018_unformatted_donor <- str_sub(martino2018_unformatted_donor,1,6)
martino2018_unformatted_donor <- as.numeric(factor(martino2018_unformatted_donor))
martino2018_formatted_donor <- paste(rep("W",205),martino2018_unformatted_donor,sep="")

martino2018_samples <- data.frame(ID=martino2018_formatted_samples,
                                  Author=rep("Martino",205),
                                  Year=rep(2018,205),
                                  Tissue=rep("Blood",205),
                                  CellType=rep("CD4+",205),
                                  Age=martino2018_formatted_ages,
                                  Condition=martino2018_formatted_condition,
                                  Sex=martino2018_formatted_sex,
                                  DonorID=martino2018_formatted_donor,
                                  Misc=rep(NA,205))

martino2018_cpgs <- read.table("Martino2018/GSE114135-GPL23976_series_matrix.txt",
                               comment = "!",
                               skip=5,
                               fill=TRUE)

martino2018_cpgs <- martino2018_cpgs[-c(1:5),]
rownames(martino2018_cpgs) <- str_sub(martino2018_cpgs[,1],1,100)
martino2018_cpgs <- martino2018_cpgs[,-1]
colnames(martino2018_cpgs) <- martino2018_formatted_samples

martino2018_cpgs <- martino2018_cpgs[rownames(martino2018_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(martino2018_cpgs),]

martino2018_cpgs$cpg <- rownames(martino2018_cpgs)
martino2018_cpgs <- data.table(martino2018_cpgs)
sample_table <- rbind(sample_table,martino2018_samples)
cpg_table <- merge(cpg_table,martino2018_cpgs, all=TRUE)

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
                                   Tissue=rep("Lung",280),
                                   CellType=rep("Airway Smooth Muscle",280),
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
thompson2020_formatted_arrangement == colnames(thompson2020_cpgs)
colnames(thompson2020_cpgs) <- thompson2020_formatted_samples

thompson2020_cpgs <- thompson2020_cpgs[rownames(thompson2020_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(thompson2020_cpgs),]

thompson2020_cpgs$cpg <- rownames(thompson2020_cpgs)
thompson2020_cpgs <- data.table(thompson2020_cpgs)
sample_table <- rbind(sample_table,thompson2020_samples)
cpg_table <- merge(cpg_table,thompson2020_cpgs, all=TRUE)


zannas2019_unformatted_table <- read.table("Zannas2019/GSE128235_series_matrix.txt",
                                           comment = "!",
                                           skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
zannas2019_unformatted_samples <- zannas2019_unformatted_table[6,-1]
zannas2019_formatted_samples <- str_sub(zannas2019_unformatted_samples,1,62)[1:537]

#Formatting ages.
zannas2019_unformatted_ages <- zannas2019_unformatted_table[3,-1]
zannas2019_formatted_ages <- as.numeric(str_sub(zannas2019_unformatted_ages, 5,100))

#Formatting sex.
zannas2019_unformatted_sex <- zannas2019_unformatted_table[4,-1]
zannas2019_formatted_sex <- str_sub(zannas2019_unformatted_sex, 6, 30)

#Formatting condition and miscellaneous information.
zannas2019_unformatted_condition <- zannas2019_unformatted_table[2,-1]
zannas2019_formatted_condition <- str_sub(zannas2019_unformatted_condition, 12, 30)
zannas2019_formatted_condition[zannas2019_formatted_condition=="control"] <- "Control"
zannas2019_formatted_condition[zannas2019_formatted_condition=="case"] <- "Depression"

zannas2019_formatted_donor <- paste(rep("Y",537),c(1:537),sep="")

zannas2019_unformatted_arrangement <- zannas2019_unformatted_table[1,-1]
zannas2019_formatted_arrangement <- str_sub(zannas2019_unformatted_arrangement,25,100)

zannas2019_samples <- data.frame(ID=zannas2019_formatted_samples,
                                 Author=rep("Zannas",537),
                                 Year=rep(2019,537),
                                 Tissue=rep("Blood",537),
                                 CellType=rep("Blood",537),
                                 Age=zannas2019_formatted_ages,
                                 Condition=zannas2019_formatted_condition,
                                 Sex=zannas2019_formatted_sex,
                                 DonorID=zannas2019_formatted_donor,
                                 Misc=rep(NA,537))

zannas2019_cpgs <- data.table::fread("Zannas2019/GSE128235_matrix_normalized.txt",
                                     header=FALSE)  %>%
  as.data.frame()
rownames(zannas2019_cpgs) <- zannas2019_cpgs[,2]
zannas2019_cpgs <- zannas2019_cpgs[,-c(1,2)]
zannas2019_cpgs <- zannas2019_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),537)]

zannas2019_colnames <- read_tsv("Zannas2019/first_line.txt",
                                col_names = FALSE)
zannas2019_colnames <- str_sub(unlist(zannas2019_colnames[1,]),1,100)
zannas2019_colnames <- zannas2019_colnames[-c(1,2)]
zannas2019_colnames <- zannas2019_colnames[rep(c(rep(TRUE, 2- 1), FALSE),537)]
zannas2019_colnames <- str_sub(zannas2019_colnames,7,15)
colnames(zannas2019_cpgs) <- zannas2019_colnames
zannas2019_cpgs <- zannas2019_cpgs[,zannas2019_formatted_arrangement]
zannas2019_formatted_arrangement == colnames(zannas2019_cpgs)
colnames(zannas2019_cpgs) <- zannas2019_formatted_samples

zannas2019_cpgs <- zannas2019_cpgs[rownames(zannas2019_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(zannas2019_cpgs),]

zannas2019_cpgs$cpg <- rownames(zannas2019_cpgs)
zannas2019_cpgs <- data.table(zannas2019_cpgs)
sample_table <- rbind(sample_table,zannas2019_samples)
cpg_table <- merge(cpg_table,zannas2019_cpgs, all=TRUE)

nicodemus2017_unformatted_table <- read.table("Nicodemus2017/GSE85566_series_matrix.txt",
                                              comment = "!",
                                              skip=0,fill=TRUE,nrows = 7)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
nicodemus2017_unformatted_samples <- nicodemus2017_unformatted_table[4,-1]
nicodemus2017_formatted_samples <- str_sub(nicodemus2017_unformatted_samples,1,62)[1:115]

#Formatting ages.
nicodemus2017_unformatted_ages <- nicodemus2017_unformatted_table[1,-1]
nicodemus2017_formatted_ages <- as.numeric(str_sub(nicodemus2017_unformatted_ages, 5,100))

#Formatting sex.
nicodemus2017_unformatted_sex <- nicodemus2017_unformatted_table[3,-1]
nicodemus2017_formatted_sex <- str_sub(nicodemus2017_unformatted_sex, 9, 30)

#Formatting condition and miscellaneous information.
nicodemus2017_unformatted_condition <- nicodemus2017_unformatted_table[2,-1]
nicodemus2017_formatted_condition <- str_sub(nicodemus2017_unformatted_condition, 17, 30)

nicodemus2017_formatted_donor <- paste(rep("Z",115),c(1:115),sep="")

nicodemus2017_samples <- data.frame(ID=nicodemus2017_formatted_samples,
                                    Author=rep("Nicodemus",115),
                                    Year=rep(2017,115),
                                    Tissue=rep("Lung",115),
                                    CellType=rep("Epithelial Cells",115),
                                    Age=nicodemus2017_formatted_ages,
                                    Condition=nicodemus2017_formatted_condition,
                                    Sex=nicodemus2017_formatted_sex,
                                    DonorID=nicodemus2017_formatted_donor,
                                    Misc=rep(NA,115))

nicodemus2017_cpgs <- read.table("Nicodemus2017/GSE85566_series_matrix.txt",
                                 comment = "!",
                                 skip=5,
                                 fill=TRUE)
nicodemus2017_cpgs <- nicodemus2017_cpgs[-c(1:4),]
rownames(nicodemus2017_cpgs) <- nicodemus2017_cpgs[,1]
nicodemus2017_cpgs <- nicodemus2017_cpgs[,-1]
colnames(nicodemus2017_cpgs) <- nicodemus2017_formatted_samples

nicodemus2017_cpgs <- nicodemus2017_cpgs[rownames(nicodemus2017_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(nicodemus2017_cpgs),]

nicodemus2017_cpgs$cpg <- rownames(nicodemus2017_cpgs)
nicodemus2017_cpgs <- data.table(nicodemus2017_cpgs)
sample_table <- rbind(sample_table,nicodemus2017_samples)
cpg_table <- merge(cpg_table,nicodemus2017_cpgs, all=TRUE)

langevin2016_unformatted_table <- read.table("Langevin2016/GSE70977_series_matrix.txt",
                                             comment = "!",
                                             skip=0,fill=TRUE,nrows = 7)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
langevin2016_unformatted_samples <- langevin2016_unformatted_table[5,-1]
langevin2016_formatted_samples <- str_sub(langevin2016_unformatted_samples,1,62)[1:223]

#Formatting ages.
langevin2016_unformatted_ages <- langevin2016_unformatted_table[2,-1]
langevin2016_formatted_ages <- as.numeric(str_sub(langevin2016_unformatted_ages, 26,100))

#Formatting sex.
langevin2016_unformatted_sex <- langevin2016_unformatted_table[3,-1]
langevin2016_formatted_sex <- str_sub(langevin2016_unformatted_sex, 6, 30)

#Formatting condition and miscellaneous information.
langevin2016_unformatted_condition <- langevin2016_unformatted_table[1,-1]
langevin2016_formatted_condition <- str_sub(langevin2016_unformatted_condition, 8, 30)

langevin2016_formatted_donor <- paste(rep("AA",223),c(1:223),sep="")

langevin2016_samples <- data.frame(ID=langevin2016_formatted_samples,
                                   Author=rep("Langevin",223),
                                   Year=rep(2016,223),
                                   Tissue=rep("Mouth",223),
                                   CellType=rep("Epithelial",223),
                                   Age=langevin2016_formatted_ages,
                                   Condition=langevin2016_formatted_condition,
                                   Sex=langevin2016_formatted_sex,
                                   DonorID=langevin2016_formatted_donor,
                                   Misc=rep(NA,223))

langevin2016_cpgs <- read.table("Langevin2016/GSE70977_series_matrix.txt",
                                comment = "!",
                                skip=7,
                                fill=TRUE)
langevin2016_cpgs <- langevin2016_cpgs[-c(1:5),]
rownames(langevin2016_cpgs) <- langevin2016_cpgs[,1]
langevin2016_cpgs <- langevin2016_cpgs[,-1]
colnames(langevin2016_cpgs) <- langevin2016_formatted_samples

langevin2016_cpgs <- langevin2016_cpgs[rownames(langevin2016_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(langevin2016_cpgs),]

langevin2016_cpgs$cpg <- rownames(langevin2016_cpgs)
langevin2016_cpgs <- data.table(langevin2016_cpgs)
sample_table <- rbind(sample_table,langevin2016_samples)
cpg_table <- merge(cpg_table,langevin2016_cpgs, all=TRUE)


pai2019_unformatted_table <- read.table("Pai2019/GSE112179_series_matrix.txt",
                                        comment = "!",
                                        skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
pai2019_unformatted_samples <- pai2019_unformatted_table[5,-1]
pai2019_formatted_samples <- str_sub(pai2019_unformatted_samples,1,62)[1:100]

#Formatting ages.
pai2019_unformatted_ages <- pai2019_unformatted_table[1,-1]
pai2019_formatted_ages <- as.numeric(str_sub(pai2019_unformatted_ages, 5,100))

#Formatting sex.
pai2019_unformatted_sex <- pai2019_unformatted_table[2,-1]
pai2019_formatted_sex <- str_sub(pai2019_unformatted_sex, 6, 30)

#Formatting condition and miscellaneous information.
pai2019_unformatted_condition <- pai2019_unformatted_table[3,-1]
pai2019_formatted_condition <- str_sub(pai2019_unformatted_condition, 10, 30)

pai2019_formatted_donor <- paste(rep("AB",100),c(1:100),sep="")

pai2019_samples <- data.frame(ID=pai2019_formatted_samples,
                              Author=rep("Pai",100),
                              Year=rep(2019,100),
                              Tissue=rep("Brain",100),
                              CellType=rep("Neurons",100),
                              Age=pai2019_formatted_ages,
                              Condition=pai2019_formatted_condition,
                              Sex=pai2019_formatted_sex,
                              DonorID=pai2019_formatted_donor,
                              Misc=rep(NA,100))

pai2019_cpgs <- read.table("Pai2019/GSE112179_series_matrix.txt",
                           comment = "!",
                           skip=5,
                           fill=TRUE)
pai2019_cpgs <- pai2019_cpgs[-c(1:5),]
rownames(pai2019_cpgs) <- pai2019_cpgs[,1]
pai2019_cpgs <- pai2019_cpgs[,-1]
colnames(pai2019_cpgs) <- pai2019_formatted_samples

pai2019_cpgs <- pai2019_cpgs[rownames(pai2019_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(pai2019_cpgs),]

pai2019_cpgs$cpg <- rownames(pai2019_cpgs)
pai2019_cpgs <- data.table(pai2019_cpgs)
sample_table <- rbind(sample_table,pai2019_samples)
cpg_table <- merge(cpg_table,pai2019_cpgs, all=TRUE)

cobben2019_unformatted_table <- read.table("Cobben2019/GSE112987_series_matrix.txt",
                                           comment = "!",
                                           skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
cobben2019_unformatted_samples <- cobben2019_unformatted_table[5,-1]
cobben2019_formatted_samples <- str_sub(cobben2019_unformatted_samples,1,62)[1:103]

#Formatting ages.
cobben2019_unformatted_ages <- cobben2019_unformatted_table[3,-1]
cobben2019_formatted_ages <- as.numeric(str_sub(cobben2019_unformatted_ages, 5,100))

#Formatting sex.
cobben2019_unformatted_sex <- cobben2019_unformatted_table[1,-1]
cobben2019_formatted_sex <- str_sub(cobben2019_unformatted_sex, 9, 30)

#Formatting condition and miscellaneous information.
cobben2019_unformatted_condition <- cobben2019_unformatted_table[2,-1]
cobben2019_formatted_condition <- str_sub(cobben2019_unformatted_condition, 16, 30)
cobben2019_formatted_condition[cobben2019_formatted_condition=="control"] <- "Control"

cobben2019_formatted_donor <- paste(rep("AC",103),c(1:103),sep="")

cobben2019_samples <- data.frame(ID=cobben2019_formatted_samples,
                                 Author=rep("Cobben",103),
                                 Year=rep(2019,103),
                                 Tissue=rep("Blood",103),
                                 CellType=rep("Blood",103),
                                 Age=cobben2019_formatted_ages,
                                 Condition=cobben2019_formatted_condition,
                                 Sex=cobben2019_formatted_sex,
                                 DonorID=cobben2019_formatted_donor,
                                 Misc=rep(NA,103))

cobben2019_cpgs <- read.table("Cobben2019/GSE112987_series_matrix.txt",
                              comment = "!",
                              skip=5,
                              fill=TRUE)
cobben2019_cpgs <- cobben2019_cpgs[-c(1:5),]
rownames(cobben2019_cpgs) <- cobben2019_cpgs[,1]
cobben2019_cpgs <- cobben2019_cpgs[,-1]
colnames(cobben2019_cpgs) <- cobben2019_formatted_samples

cobben2019_cpgs <- cobben2019_cpgs[rownames(cobben2019_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(cobben2019_cpgs),]

cobben2019_cpgs$cpg <- rownames(cobben2019_cpgs)
cobben2019_cpgs <- data.table(cobben2019_cpgs)
sample_table <- rbind(sample_table,cobben2019_samples)
cpg_table <- merge(cpg_table,cobben2019_cpgs, all=TRUE)


husquin2019_unformatted_table <- read.table("Husquin2019/GSE120610_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
husquin2019_unformatted_samples <- husquin2019_unformatted_table[5,-1]
husquin2019_formatted_samples <- str_sub(husquin2019_unformatted_samples,1,62)[1:156]

#Formatting ages.
husquin2019_unformatted_ages <- husquin2019_unformatted_table[1,-1]
husquin2019_formatted_ages <- as.numeric(str_sub(husquin2019_unformatted_ages, 5,100))

#Formatting sex.
husquin2019_unformatted_sex <- husquin2019_unformatted_table[2,-1]
husquin2019_formatted_sex <- str_sub(husquin2019_unformatted_sex, 9, 15)

#Formatting condition and miscellaneous information.
husquin2019_formatted_condition <- rep("Control",156)

husquin2019_formatted_donor <- paste(rep("AD",156),c(1:156),sep="")

husquin2019_unformatted_arrangement <- husquin2019_unformatted_table[3,-1]
husquin2019_formatted_arrangement <- str_sub(husquin2019_unformatted_arrangement,8,100)

husquin2019_samples <- data.frame(ID=husquin2019_formatted_samples,
                                  Author=rep("Husquin",156),
                                  Year=rep(2019,156),
                                  Tissue=rep("Blood",156),
                                  CellType=rep("Monocytes",156),
                                  Age=husquin2019_formatted_ages,
                                  Condition=husquin2019_formatted_condition,
                                  Sex=husquin2019_formatted_sex,
                                  DonorID=husquin2019_formatted_donor,
                                  Misc=rep(NA,156))

husquin2019_cpgs <- data.table::fread("Husquin2019/GSE120610_Matrix_processed.csv",
                                      header=FALSE)  %>%
  as.data.frame()
row.names(husquin2019_cpgs) <- husquin2019_cpgs$V1
husquin2019_cpgs <- husquin2019_cpgs[,-1]
husquin2019_cpgs <- husquin2019_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),156)]

colnames(husquin2019_cpgs) <- str_sub(husquin2019_cpgs[1,],7,15)
colnames(husquin2019_cpgs) <- husquin2019_formatted_samples
husquin2019_cpgs <- husquin2019_cpgs[-1,]

husquin2019_cpgs <- husquin2019_cpgs[rownames(husquin2019_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(husquin2019_cpgs),]

husquin2019_cpgs$cpg <- rownames(husquin2019_cpgs)
husquin2019_cpgs <- data.table(husquin2019_cpgs)
sample_table <- rbind(sample_table,husquin2019_samples)
cpg_table <- merge(cpg_table,husquin2019_cpgs, all=TRUE)

gasparoni2018_unformatted_table <- read.table("Gasparoni2018/GSE66351_series_matrix.txt",
                                              comment = "!",
                                              skip=0,fill=TRUE,nrows = 8)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
gasparoni2018_unformatted_samples <- gasparoni2018_unformatted_table[8,-1]
gasparoni2018_formatted_samples <- str_sub(gasparoni2018_unformatted_samples,1,62)[1:190]

#Formatting ages.
gasparoni2018_unformatted_ages <- gasparoni2018_unformatted_table[4,-1]
gasparoni2018_formatted_ages <- as.numeric(str_sub(gasparoni2018_unformatted_ages, 5,100))

#Formatting sex.
gasparoni2018_unformatted_sex <- gasparoni2018_unformatted_table[5,-1]
gasparoni2018_formatted_sex <- str_sub(gasparoni2018_unformatted_sex, 6, 15)

#Formatting condition and miscellaneous information.
gasparoni2018_unformatted_condition <- gasparoni2018_unformatted_table[2,-1]
gasparoni2018_unformatted_condition <- str_sub(gasparoni2018_unformatted_condition,12,20)
gasparoni2018_unformatted_condition[gasparoni2018_unformatted_condition=="CTRL"] <- "Control"
gasparoni2018_unformatted_condition[gasparoni2018_unformatted_condition=="AD"] <- "Alzheimers"
gasparoni2018_formatted_condition <- gasparoni2018_unformatted_condition

gasparoni2018_unformatted_donor <- gasparoni2018_unformatted_table[6,-1]
gasparoni2018_unformatted_donor <- str_sub(gasparoni2018_unformatted_donor,11,30)
gasparoni2018_unformatted_donor <- as.numeric(factor(gasparoni2018_unformatted_donor))
gasparoni2018_formatted_donor <- paste(rep("AE",190),gasparoni2018_unformatted_donor,sep="")

gasparoni2018_unformatted_misc <- gasparoni2018_unformatted_table[3,-1]
gasparoni2018_formatted_misc <- str_sub(gasparoni2018_unformatted_misc,15,200)

gasparoni2018_samples <- data.frame(ID=gasparoni2018_formatted_samples,
                                    Author=rep("Gasparoni",190),
                                    Year=rep(2018,190),
                                    Tissue=rep("Brain",190),
                                    CellType=rep("Brain",190),
                                    Age=gasparoni2018_formatted_ages,
                                    Condition=gasparoni2018_formatted_condition,
                                    Sex=gasparoni2018_formatted_sex,
                                    DonorID=gasparoni2018_formatted_donor,
                                    Misc=gasparoni2018_formatted_misc)

gasparoni2018_cpgs <- read.table("Gasparoni2018/GSE66351_series_matrix.txt",
                                 comment = "!",
                                 skip=8,
                                 fill=TRUE)
gasparoni2018_cpgs <- gasparoni2018_cpgs[-c(1:8),]
rownames(gasparoni2018_cpgs) <- gasparoni2018_cpgs[,1]
gasparoni2018_cpgs <- gasparoni2018_cpgs[,-1]
colnames(gasparoni2018_cpgs) <- gasparoni2018_formatted_samples

gasparoni2018_cpgs <- gasparoni2018_cpgs[rownames(gasparoni2018_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(gasparoni2018_cpgs),]

gasparoni2018_cpgs$cpg <- rownames(gasparoni2018_cpgs)
gasparoni2018_cpgs <- data.table(gasparoni2018_cpgs)
sample_table <- rbind(sample_table,gasparoni2018_samples)
cpg_table <- merge(cpg_table,gasparoni2018_cpgs, all=TRUE)


liu2018_unformatted_table <- read.table("Liu2018/GSE106648_series_matrix.txt",
                                        comment = "!",
                                        skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
liu2018_unformatted_samples <- liu2018_unformatted_table[5,-1]
liu2018_formatted_samples <- str_sub(liu2018_unformatted_samples,1,62)[1:279]

#Formatting ages.
liu2018_unformatted_ages <- liu2018_unformatted_table[3,-1]
liu2018_formatted_ages <- as.numeric(str_sub(liu2018_unformatted_ages, 5,100))

#Formatting sex.
liu2018_unformatted_sex <- liu2018_unformatted_table[2,-1]
liu2018_formatted_sex <- str_sub(liu2018_unformatted_sex, 9, 15)

#Formatting condition and miscellaneous information.
liu2018_unformatted_condition <- liu2018_unformatted_table[1,-1]
liu2018_unformatted_condition <- str_sub(liu2018_unformatted_condition,17,40)
liu2018_unformatted_condition[liu2018_unformatted_condition=="Healthy control"] = "Control"
liu2018_unformatted_condition[liu2018_unformatted_condition=="MS case"] = "MS"
liu2018_formatted_condition <- liu2018_unformatted_condition

liu2018_formatted_donor <- paste(rep("AF",279),c(279),sep="")

liu2018_samples <- data.frame(ID=liu2018_formatted_samples,
                              Author=rep("Liu",279),
                              Year=rep(2018,279),
                              Tissue=rep("Blood",279),
                              CellType=rep("Blood",279),
                              Age=liu2018_formatted_ages,
                              Condition=liu2018_formatted_condition,
                              Sex=liu2018_formatted_sex,
                              DonorID=liu2018_formatted_donor,
                              Misc=rep(NA,279))

liu2018_cpgs <- read.table("Liu2018/GSE106648_series_matrix.txt",
                           comment = "!",
                           skip=5,
                           fill=TRUE)
liu2018_cpgs <- liu2018_cpgs[-c(1:5),]
rownames(liu2018_cpgs) <- liu2018_cpgs[,1]
liu2018_cpgs <- liu2018_cpgs[,-1]
colnames(liu2018_cpgs) <- liu2018_formatted_samples

liu2018_cpgs <- liu2018_cpgs[rownames(liu2018_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(liu2018_cpgs),]

liu2018_cpgs$cpg <- rownames(liu2018_cpgs)
liu2018_cpgs <- data.table(liu2018_cpgs)
sample_table <- rbind(sample_table,liu2018_samples)
cpg_table <- merge(cpg_table,liu2018_cpgs, all=TRUE)


somineni2018_unformatted_table <- read.table("Somineni2018/GSE112611_series_matrix.txt",
                                             comment = "!",
                                             skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
somineni2018_unformatted_samples <- somineni2018_unformatted_table[6,-1]
somineni2018_formatted_samples <- str_sub(somineni2018_unformatted_samples,1,62)[1:402]

#Formatting ages.
somineni2018_unformatted_ages <- somineni2018_unformatted_table[2,-1]
somineni2018_formatted_ages <- as.numeric(str_sub(somineni2018_unformatted_ages, 5,100))

#Formatting sex.
somineni2018_unformatted_sex <- somineni2018_unformatted_table[1,-1]
somineni2018_formatted_sex <- str_sub(somineni2018_unformatted_sex, 9, 15)

#Formatting condition and miscellaneous information.
somineni2018_unformatted_condition <- somineni2018_unformatted_table[3,-1]
somineni2018_unformatted_condition <- str_sub(somineni2018_unformatted_condition,12,40)
somineni2018_unformatted_condition[somineni2018_unformatted_condition=="non-IBD control"] = "Control"
somineni2018_formatted_condition <- somineni2018_unformatted_condition

somineni2018_formatted_donor <- paste(rep("AG",402),c(1:402),sep="")

somineni2018_formatted_location <- str_sub(somineni2018_unformatted_table[4,-1],1,100)

somineni2018_samples <- data.frame(ID=somineni2018_formatted_samples,
                                   Author=rep("Somineni",402),
                                   Year=rep(2018,402),
                                   Tissue=rep("Blood",402),
                                   CellType=rep("Blood",402),
                                   Age=somineni2018_formatted_ages,
                                   Condition=somineni2018_formatted_condition,
                                   Sex=somineni2018_formatted_sex,
                                   DonorID=somineni2018_formatted_donor,
                                   Misc=rep(NA,402))

somineni2018_cpgs <- data.table::fread("Somineni2018/GSE112611_beta_values.txt",
                                       header=FALSE)  %>%
  as.data.frame()
row.names(somineni2018_cpgs) <- somineni2018_cpgs$V1
somineni2018_cpgs <- somineni2018_cpgs[,-1]
somineni2018_cpgs <- somineni2018_cpgs[, rep(c(rep(TRUE, 2- 1), FALSE),402)]
colnames(somineni2018_cpgs) <- somineni2018_cpgs[1,]
somineni2018_cpgs <- somineni2018_cpgs[-1,]

colnames(somineni2018_cpgs) <- somineni2018_formatted_samples

somineni2018_cpgs <- somineni2018_cpgs[rownames(somineni2018_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(somineni2018_cpgs),]

somineni2018_cpgs$cpg <- rownames(somineni2018_cpgs)
somineni2018_cpgs <- data.table(somineni2018_cpgs)
sample_table <- rbind(sample_table,somineni2018_samples)
cpg_table <- merge(cpg_table,somineni2018_cpgs, all=TRUE)

roos2017_unformatted_table <- read.table("Roos2017/GSE90124_series_matrix.txt",
                                         comment = "!",
                                         skip=0,fill=TRUE,nrows = 6)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
roos2017_unformatted_samples <- roos2017_unformatted_table[3,-1]
roos2017_formatted_samples <- str_sub(roos2017_unformatted_samples,1,62)[1:322]

#Formatting ages.
roos2017_unformatted_ages <- roos2017_unformatted_table[1,-1]
roos2017_formatted_ages <- as.numeric(str_sub(roos2017_unformatted_ages, 15,100))

#Formatting sex.
roos2017_formatted_sex <- rep("Female",322)

#Formatting condition and miscellaneous information.
roos2017_formatted_condition <- rep("Control",322)

roos2017_formatted_donor <- paste(rep("AH",322),c(1:322),sep="")

roos2017_samples <- data.frame(ID=roos2017_formatted_samples,
                               Author=rep("Roos",322),
                               Year=rep(2017,322),
                               Tissue=rep("Skin",322),
                               CellType=rep("Epithelial",322),
                               Age=roos2017_formatted_ages,
                               Condition=roos2017_formatted_condition,
                               Sex=roos2017_formatted_sex,
                               DonorID=roos2017_formatted_donor,
                               Misc=rep(NA,322))

roos2017_cpgs <- read.table("Roos2017/GSE90124_series_matrix.txt",
                            comment = "!",
                            skip=3,
                            fill=TRUE)
roos2017_cpgs <- roos2017_cpgs[-c(1:3),]
rownames(roos2017_cpgs) <- roos2017_cpgs[,1]
roos2017_cpgs <- roos2017_cpgs[,-1]
colnames(roos2017_cpgs) <- roos2017_formatted_samples

roos2017_cpgs <- roos2017_cpgs[rownames(roos2017_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(roos2017_cpgs),]

roos2017_cpgs$cpg <- rownames(roos2017_cpgs)
roos2017_cpgs <- data.table(roos2017_cpgs)
sample_table <- rbind(sample_table,roos2017_samples)
cpg_table <- merge(cpg_table,roos2017_cpgs, all=TRUE)


kananen2016_unformatted_table <- read.table("Kananen2016/GSE69270_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 4)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
kananen2016_unformatted_samples <- kananen2016_unformatted_table[4,-1]
kananen2016_formatted_samples <- str_sub(kananen2016_unformatted_samples,1,62)[1:184]

#Formatting ages.
kananen2016_unformatted_ages <- kananen2016_unformatted_table[1,-1]
kananen2016_formatted_ages <- as.numeric(str_sub(kananen2016_unformatted_ages, 14,100))

#Formatting sex.
kananen2016_unformatted_sex <- kananen2016_unformatted_table[2,-1]
kananen2016_unformatted_sex <- str_sub(kananen2016_unformatted_sex,-1,-1)
kananen2016_unformatted_sex[kananen2016_unformatted_sex=="1"] <- "Female"
kananen2016_unformatted_sex[kananen2016_unformatted_sex=="0"] <- "Male"
kananen2016_formatted_sex <- kananen2016_unformatted_sex

#Formatting condition and miscellaneous information.
kananen2016_formatted_condition <- rep("Control",184)

kananen2016_formatted_donor <- paste(rep("AI",184),c(1:184),sep="")

kananen2016_samples <- data.frame(ID=kananen2016_formatted_samples,
                                  Author=rep("Kananen",184),
                                  Year=rep(2016,184),
                                  Tissue=rep("Blood",184),
                                  CellType=rep("Leukocytes",184),
                                  Age=kananen2016_formatted_ages,
                                  Condition=kananen2016_formatted_condition,
                                  Sex=kananen2016_formatted_sex,
                                  DonorID=kananen2016_formatted_donor,
                                  Misc=rep(NA,184))

kananen2016_cpgs <- read.table("Kananen2016/GSE69270_series_matrix.txt",
                               comment = "!",
                               skip=4,
                               fill=TRUE)
kananen2016_cpgs <- kananen2016_cpgs[-c(1:4),]
rownames(kananen2016_cpgs) <- kananen2016_cpgs[,1]
kananen2016_cpgs <- kananen2016_cpgs[,-1]
colnames(kananen2016_cpgs) <- kananen2016_formatted_samples

kananen2016_cpgs <- kananen2016_cpgs[rownames(kananen2016_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(kananen2016_cpgs),]

kananen2016_cpgs$cpg <- rownames(kananen2016_cpgs)
kananen2016_cpgs <- data.table(kananen2016_cpgs)
sample_table <- rbind(sample_table,kananen2016_samples)
cpg_table <- merge(cpg_table,kananen2016_cpgs, all=TRUE)

cerapio2021_unformatted_table <- read.table("Cerapio2021/GSE136583_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
cerapio2021_unformatted_samples <- cerapio2021_unformatted_table[5,-1]
cerapio2021_formatted_samples <- str_sub(cerapio2021_unformatted_samples,1,62)[1:62]

#Formatting ages.
cerapio2021_unformatted_ages <- cerapio2021_unformatted_table[2,-1]
cerapio2021_formatted_ages <- as.numeric(str_sub(cerapio2021_unformatted_ages, 5,100))

#Formatting sex.
cerapio2021_unformatted_sex <- cerapio2021_unformatted_table[3,-1]
cerapio2021_unformatted_sex <- str_sub(cerapio2021_unformatted_sex,9,30)
cerapio2021_formatted_sex <- cerapio2021_unformatted_sex

#Formatting condition and miscellaneous information.
cerapio2021_unformatted_condition <- cerapio2021_unformatted_table[1,-1]
cerapio2021_unformatted_condition <- str_sub(cerapio2021_unformatted_condition,9,50)
cerapio2021_unformatted_condition[cerapio2021_unformatted_condition=="non-tumor liver"] <- "Control"
cerapio2021_formatted_condition <- cerapio2021_unformatted_condition

cerapio2021_formatted_donor <- paste(rep("AJ",62),c(1:62),sep="")

cerapio2021_samples <- data.frame(ID=cerapio2021_formatted_samples,
                                  Author=rep("Cerapio",62),
                                  Year=rep(2021,62),
                                  Tissue=rep("Liver",62),
                                  CellType=rep("Hepatocytes",62),
                                  Age=cerapio2021_formatted_ages,
                                  Condition=cerapio2021_formatted_condition,
                                  Sex=cerapio2021_formatted_sex,
                                  DonorID=cerapio2021_formatted_donor,
                                  Misc=rep(NA,62))

cerapio2021_cpgs <- read.table("Cerapio2021/GSE136583_series_matrix.txt",
                               comment = "!",
                               skip=5,
                               fill=TRUE)
cerapio2021_cpgs <- cerapio2021_cpgs[-c(1:5),]
rownames(cerapio2021_cpgs) <- cerapio2021_cpgs[,1]
cerapio2021_cpgs <- cerapio2021_cpgs[,-1]
colnames(cerapio2021_cpgs) <- cerapio2021_formatted_samples

cerapio2021_cpgs <- cerapio2021_cpgs[rownames(cerapio2021_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(cerapio2021_cpgs),]

cerapio2021_cpgs$cpg <- rownames(cerapio2021_cpgs)
cerapio2021_cpgs <- data.table(cerapio2021_cpgs)
sample_table <- rbind(sample_table,cerapio2021_samples)
cpg_table <- merge(cpg_table,cerapio2021_cpgs, all=TRUE)

hearn2020_unformatted_table <- read.table("Hearn2020/GSE119078_series_matrix.txt",
                                          comment = "!",
                                          skip=0,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
hearn2020_unformatted_samples <- hearn2020_unformatted_table[5,-1]
hearn2020_formatted_samples <- str_sub(hearn2020_unformatted_samples,1,62)[1:59]

#Formatting ages.
hearn2020_unformatted_ages <- hearn2020_unformatted_table[2,-1]
hearn2020_formatted_ages <- as.numeric(str_sub(hearn2020_unformatted_ages, 5,100))

#Formatting sex.
hearn2020_unformatted_sex <- hearn2020_unformatted_table[1,-1]
hearn2020_unformatted_sex <- str_sub(hearn2020_unformatted_sex,9,100)
hearn2020_formatted_sex <- hearn2020_unformatted_sex

#Formatting condition and miscellaneous information.
hearn2020_unformatted_condition <- hearn2020_unformatted_table[3,-1]
hearn2020_unformatted_condition <- str_sub(hearn2020_unformatted_condition,17,40)
hearn2020_formatted_condition <- hearn2020_unformatted_condition

hearn2020_formatted_donor <- paste(rep("AK",59),c(1:59),sep="")

hearn2020_samples <- data.frame(ID=hearn2020_formatted_samples,
                                Author=rep("Hearn",59),
                                Year=rep(2020,59),
                                Tissue=rep("Saliva",59),
                                CellType=rep("Epithelial",59),
                                Age=hearn2020_formatted_ages,
                                Condition=hearn2020_formatted_condition,
                                Sex=hearn2020_formatted_sex,
                                DonorID=hearn2020_formatted_donor,
                                Misc=rep(NA,59))

hearn2020_cpgs <- read.table("Hearn2020/GSE119078_series_matrix.txt",
                             comment = "!",
                             skip=5,
                             fill=TRUE)
hearn2020_cpgs <- hearn2020_cpgs[-c(1:4),]
rownames(hearn2020_cpgs) <- hearn2020_cpgs[,1]
hearn2020_cpgs <- hearn2020_cpgs[,-1]
colnames(hearn2020_cpgs) <- hearn2020_formatted_samples

hearn2020_cpgs <- hearn2020_cpgs[rownames(hearn2020_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(hearn2020_cpgs),]

hearn2020_cpgs$cpg <- rownames(hearn2020_cpgs)
hearn2020_cpgs <- data.table(hearn2020_cpgs)
sample_table <- rbind(sample_table,hearn2020_samples)
cpg_table <- merge(cpg_table,hearn2020_cpgs, all=TRUE)

hong2017_unformatted_table <- read.table("Hong2017/GSE92767_series_matrix.txt",
                                         comment = "!",
                                         skip=0,fill=TRUE,nrows = 4)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
hong2017_unformatted_samples <- hong2017_unformatted_table[4,-1]
hong2017_formatted_samples <- str_sub(hong2017_unformatted_samples,1,62)[1:54]

#Formatting ages.
hong2017_unformatted_ages <- hong2017_unformatted_table[1,-1]
hong2017_formatted_ages <- as.numeric(str_sub(hong2017_unformatted_ages, 6,100))

#Formatting sex.
hong2017_unformatted_sex <- hong2017_unformatted_table[2,-1]
hong2017_unformatted_sex <- str_sub(hong2017_unformatted_sex,9,15)
hong2017_formatted_sex <- hong2017_unformatted_sex

#Formatting condition and miscellaneous information.
hong2017_formatted_condition <- rep("Control",54)

hong2017_formatted_donor <- paste(rep("AL",54),c(1:54),sep="")

hong2017_samples <- data.frame(ID=hong2017_formatted_samples,
                               Author=rep("Hong",54),
                               Year=rep(2017,54),
                               Tissue=rep("Saliva",54),
                               CellType=rep("Epithelial",54),
                               Age=hong2017_formatted_ages,
                               Condition=hong2017_formatted_condition,
                               Sex=hong2017_formatted_sex,
                               DonorID=hong2017_formatted_donor,
                               Misc=rep(NA,54))

hong2017_cpgs <- read.table("Hong2017/GSE92767_series_matrix.txt",
                            comment = "!",
                            skip=4,
                            fill=TRUE)
hong2017_cpgs <- hong2017_cpgs[-c(1:4),]
rownames(hong2017_cpgs) <- hong2017_cpgs[,1]
hong2017_cpgs <- hong2017_cpgs[,-1]
colnames(hong2017_cpgs) <- hong2017_formatted_samples

hong2017_cpgs <- hong2017_cpgs[rownames(hong2017_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(hong2017_cpgs),]

hong2017_cpgs$cpg <- rownames(hong2017_cpgs)
hong2017_cpgs <- data.table(hong2017_cpgs)
sample_table <- rbind(sample_table,hong2017_samples)
cpg_table <- merge(cpg_table,hong2017_cpgs, all=TRUE)


gopalan2017_unformatted_table <- read.table("Gopalan2017/GSE99091-GPL13534_series_matrix.txt",
                                            comment = "!",
                                            skip=0,fill=TRUE,nrows = 4)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
gopalan2017_unformatted_samples <- gopalan2017_unformatted_table[4,-1]
gopalan2017_formatted_samples <- str_sub(gopalan2017_unformatted_samples,1,62)[1:57]

#Formatting ages.
gopalan2017_unformatted_ages <- gopalan2017_unformatted_table[2,-1]
gopalan2017_formatted_ages <- as.numeric(str_sub(gopalan2017_unformatted_ages, 6,100))

#Formatting sex.
gopalan2017_unformatted_sex <- gopalan2017_unformatted_table[1,-1]
gopalan2017_unformatted_sex <- str_sub(gopalan2017_unformatted_sex,9,60)
gopalan2017_unformatted_sex[gopalan2017_unformatted_sex=="female"] <- "Female"
gopalan2017_unformatted_sex[gopalan2017_unformatted_sex=="male"] <- "Male"
gopalan2017_formatted_sex <- gopalan2017_unformatted_sex

#Formatting condition and miscellaneous information.
gopalan2017_formatted_condition <- rep("Control",57)

gopalan2017_formatted_donor <- paste(rep("AM",57),c(1:57),sep="")

gopalan2017_samples <- data.frame(ID=gopalan2017_formatted_samples,
                                  Author=rep("Gopalan",57),
                                  Year=rep(2017,57),
                                  Tissue=rep("Saliva",57),
                                  CellType=rep("Epithelial",57),
                                  Age=gopalan2017_formatted_ages,
                                  Condition=gopalan2017_formatted_condition,
                                  Sex=gopalan2017_formatted_sex,
                                  DonorID=gopalan2017_formatted_donor,
                                  Misc=rep(NA,57))

gopalan2017_cpgs <- read.table("Gopalan2017/GSE99091-GPL13534_series_matrix.txt",
                               comment = "!",
                               skip=4,
                               fill=TRUE)
gopalan2017_cpgs <- gopalan2017_cpgs[-c(1:4),]
rownames(gopalan2017_cpgs) <- gopalan2017_cpgs[,1]
gopalan2017_cpgs <- gopalan2017_cpgs[,-1]
colnames(gopalan2017_cpgs) <- gopalan2017_formatted_samples

gopalan2017_cpgs <- gopalan2017_cpgs[rownames(gopalan2017_cpgs) %in% shared_cpgs,]
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(gopalan2017_cpgs),]

gopalan2017_cpgs$cpg <- rownames(gopalan2017_cpgs)
gopalan2017_cpgs <- data.table(gopalan2017_cpgs)
sample_table <- rbind(sample_table,gopalan2017_samples)
cpg_table <- merge(cpg_table,gopalan2017_cpgs, all=TRUE)

###########################################################################

cpg_table <- na.omit(cpg_table)
sample_table <- sample_table[sample_table$ID %in% colnames(cpg_table),]
cpg_table <- cpg_table[,colnames(cpg_table) %in% sample_table$ID]
cpg_table <- cpg_table[,unique(colnames(cpg_table))]

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
