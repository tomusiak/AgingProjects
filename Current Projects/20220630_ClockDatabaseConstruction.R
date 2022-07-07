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

sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
cpg_table <- read.csv("ClockConstruction/cpg_table.csv",row.names=1)

#Reading in the 450K dataset from Magnaye 2022 et al.
magnaye2022450k_unformatted_table <- read.table("Magnaye2022/GSE201872-GPL21145_series_matrix.txt",
                                            comment = "!",
                                            skip=5,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
magnaye2022450k_unformatted_samples <- magnaye2022450k_unformatted_table[2,-1]
magnaye2022450k_formatted_samples <- t(magnaye2022450k_unformatted_samples)[,1]

#Formatting ages.
magnaye2022450k_unformatted_ages <- magnaye2022450k_unformatted_table[3,-1]
magnaye2022450k_formatted_ages <- as.numeric(str_sub(magnaye2022450k_unformatted_ages, 6, 7))

#Formatting sex.
magnaye2022450k_unformatted_sex <- magnaye2022450k_unformatted_table[4,-1]
magnaye2022450k_formatted_sex <- str_sub(magnaye2022450k_unformatted_sex, 6, 15)

#Formatting condition.
magnaye2022450k_unformatted_condition <- magnaye2022450k_unformatted_table[5,-1]
magnaye2022450k_formatted_condition <- str_sub(magnaye2022450k_unformatted_condition, 9, 18)

#Looks like the values provided are heavily pre-processed m-values. To ensure consistency,
# I will process raw data to generate normalized beta values.
magnaye2022450k_list_of_files<-list.files(file.path("Magnaye2022/raw_data"))
magnaye2022450k_list_of_files <- magnaye2022450k_list_of_files[231:length(magnaye2022450k_list_of_files)]
magnaye2022450k_parsed_list_of_files <- substr(magnaye2022450k_list_of_files,1,30)
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

#Reading in the EPIC dataset from Magnaye 2022 et al.
magnaye2022EPIC_unformatted_table <- read.table("Magnaye2022/GSE201872-GPL13534_series_matrix.txt",
                                                comment = "!",
                                                skip=5,fill=TRUE,nrows = 5)

#Formatting sample names, as the header text file is not arranged in a particularly easy way.
magnaye2022EPIC_unformatted_samples <- magnaye2022EPIC_unformatted_table[1,-1]
magnaye2022EPIC_formatted_samples <- t(magnaye2022EPIC_unformatted_samples)[,1]

#Formatting ages.
magnaye2022EPIC_unformatted_ages <- magnaye2022EPIC_unformatted_table[2,-1]
magnaye2022EPIC_formatted_ages <- as.numeric(str_sub(magnaye2022EPIC_unformatted_ages, 6, 7))

#Formatting sex.
magnaye2022EPIC_unformatted_sex <- magnaye2022EPIC_unformatted_table[3,-1]
magnaye2022EPIC_formatted_sex <- str_sub(magnaye2022EPIC_unformatted_sex, 6, 15)

#Formatting condition.
magnaye2022EPIC_unformatted_condition <- magnaye2022EPIC_unformatted_table[4,-1]
magnaye2022EPIC_formatted_condition <- str_sub(magnaye2022EPIC_unformatted_condition, 9, 18)

#Looks like the values provided are heavily pre-processed m-values. To ensure consistency,
# I will process raw data to generate normalized beta values.
magnaye2022EPIC_list_of_files<-list.files(file.path("Magnaye2022/raw_data"))
magnaye2022EPIC_list_of_files <- magnaye2022EPIC_list_of_files[7:230]
magnaye2022EPIC_parsed_list_of_files <- substr(magnaye2022EPIC_list_of_files,1,28)
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
magnaye2022EPIC_cpgs <- magnaye2022EPIC_cpgs[rownames(magnaye2022EPIC_cpgs) %in%
                                               rownames(magnaye2022450k_cpgs),]
magnaye2022450k_cpgs <- magnaye2022450k_cpgs[rownames(magnaye2022450k_cpgs) %in%
                                               rownames(magnaye2022EPIC_cpgs),]

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
magnayefirst_samples <- data.frame(ID=magnaye2022450k_formatted_samples,
                                   Author=rep("Magnaye",30),
                                   Year=rep(2022,30),
                                   Tissue=rep("Bronchi",30),
                                   CellType=rep("Epithelial",30),
                                   Age=magnaye2022450k_formatted_ages,
                                   Condition=magnaye2022450k_formatted_condition,
                                   Sex=magnaye2022450k_formatted_sex,
                                   DonorID=magnayefirst_IDs,
                                   Misc=rep("",30))
magnayesecond_IDs <- paste(rep("G",112),1:112,sep="")
magnayesecond_samples <- data.frame(ID=magnaye2022EPIC_formatted_samples,
                                   Author=rep("Magnaye",112),
                                   Year=rep(2022,112),
                                   Tissue=rep("Bronchi",112),
                                   CellType=rep("Epithelial",112),
                                   Age=magnaye2022EPIC_formatted_ages,
                                   Condition=magnaye2022EPIC_formatted_condition,
                                   Sex=magnaye2022EPIC_formatted_sex,
                                   DonorID=magnayesecond_IDs,
                                   Misc=rep("",112))
sample_table <- rbind(magnayefirst_samples,magnayesecond_samples)
cpg_table <- cbind(magnaye2022450k_cpgs,magnaye2022EPIC_cpgs)


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
estupinan2022_cpgs <- estupinan2022_cpgs[rownames(estupinan2022_cpgs) %in% rownames(cpg_table),]
cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(estupinan2022_cpgs),]

sample_table <- rbind(sample_table,estupinan2022_samples)
cpg_table <- cbind(cpg_table,estupinan2022_cpgs)


#Time for a blood data set.
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
okereke2021_cpgs <- okereke2021_cpgs[-c(1:5)]
rownames(okereke2021_cpgs) <- okereke2021_cpgs$V1
okereke2021_cpgs <- data.frame(okereke2021_cpgs[,-1])
colnames(okereke2021_cpgs) <- okereke2021_unformatted_samples
okereke2021_cpgs <- okereke2021_cpgs[rownames(okereke2021_cpgs) %in% rownames(cpg_table),]
cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(okereke2021_cpgs),]

sample_table <- rbind(sample_table,okereke2021_samples)
cpg_table <- cbind(cpg_table,okereke2021_cpgs)

# Skin data set.
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
muse2021_cpgs <- muse2021_cpgs[rownames(muse2021_cpgs) %in% rownames(cpg_table),]
cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(muse2021_cpgs),]

sample_table <- rbind(sample_table,muse2021_samples)
cpg_table <- cbind(cpg_table,muse2021_cpgs)
write.csv(cpg_table,"ClockConstruction/cpg_table.csv")
write.csv(sample_table,"ClockConstruction/sample_table.csv")
