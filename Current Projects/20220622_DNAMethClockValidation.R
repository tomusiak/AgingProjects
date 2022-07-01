#In this program I will construct a database of DNA methylation data, consisting of 450K
# chip results. I will construct this database as having two segments - "cpg_table" consisting
# of cpgs in rows and samples in columns, and "sample_table" consisting of samples in rows and
# metadata in columns.

source("AgingProjects/Useful Scripts/generally_useful.R") #Helper functions

#Packages
setwd("Data/") #Sets directory.
library(readr)
library(tidyr)
library(dplyr)
library(devtools)
library(ggplot2)
library("methylclock")

#Quick helper function for pre-processing datasets.
# preprocessDataset <- function(file_name) {
#   cpgs <- data.frame(read_table2(file_name,
#                                  col_names = FALSE, comment = "!"))
#   colnames(cpgs) <- cpgs[1,]
#   rownames(cpgs) <- cpgs[,1]
#   cpgs <- data.frame(cpgs[-1,-1])
#   cpgs <- mutate_all(cpgs, function(x) as.numeric(as.character(x)))
#   cpgs <- drop_na(cpgs)
#   rownames(cpgs) <- str_replace_all(rownames(cpgs), "[[:punct:]]", " ")
#   rownames(cpgs) <- gsub(" ","",rownames(cpgs))
#   colnames(cpgs) <- substr(colnames(cpgs),3,12)
#   return(cpgs)
# }
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
# 
# #Pre-processing the Garaud 2017 by structuring it appropriately and annotating
# # each sample.
# garaud_cpgs <- preprocessDataset("Garaud2017/GSE71825_series_matrix.txt")
# garaud_samples <- data.frame(ID=colnames(garaud_cpgs),
#                              Author=rep("Garaud",12),
#                              Year=rep(2017,12),
#                              Tissue=rep("Blood",12),
#                              CellType=c(rep("Naive CD4+",3),rep("Memory CD4+",3),
#                                         rep("Naive CD4+",3),rep("Memory CD4+",3)),
#                              Age=rep(NA,12),
#                              Condition=rep("Healthy",12),
#                              Sex=rep(NA,12),
#                              DonorID=c("A1","A2","A3","A1","A2","A3",
#                                        "A1","A2","A3","A1","A2","A3"),
#                              Misc=c(rep("Resting",6),rep("Activated",6)))
# 
# #Initializing  the cpg_table with the first dataset.
# cpg_table <- garaud_cpgs
# 
# #Adding Garaud 2017 to the data table.
# sample_table <- bind_rows(sample_table,garaud_samples)
# 
# #Processing and adding Schlums 2015
# schlums_cpgs <- preprocessDataset("Schlums2015/GSE66564-GPL13534_series_matrix.txt")
# schlums_cpgs <- schlums_cpgs[,c(2,5,7,9,13,16,19,23)] #Sorting out NK cells.
# schlums_samples <- data.frame(ID=colnames(schlums_cpgs),
#                               Author=rep("Schlums",8),
#                               Year=rep(2015,8),
#                               Tissue=rep("Blood",8),
#                               CellType=c(rep(c("Effector CD8+","Naive CD8+"),4)),
#                               Age=rep(NA,8),
#                               Condition=rep("Healthy",8),
#                               Sex=rep(NA,8),
#                               DonorID=c("B1","B1","B2","B2",
#                                         "B3","B3","B4","B4"),
#                               Misc=c(rep("Resting",8)))
# sample_table <- rbind(sample_table,schlums_samples)
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(schlums_cpgs),]
# schlums_cpgs <- schlums_cpgs[rownames(schlums_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cbind(cpg_table,schlums_cpgs)
# 
# #Adding Rodriguez 2017 to the sample table.
# rodriguez_cpgs <- preprocessDataset("Rodriguez2017/GSE83159-GPL13534_series_matrix.txt")
# rodriguez_samples <- data.frame(ID=colnames(rodriguez_cpgs),
#                                 Author=rep("Rodriguez",6),
#                                 Year=rep(2017,6),
#                                 Tissue=rep("Blood",6),
#                                 CellType=c("Naive CD8+","Naive CD8+",
#                                            "TEMRA CD8+","TEMRA CD8+",
#                                            "Effector CD8+","Effector CD8+"),
#                                 Age=rep(NA,6),
#                                 Condition=rep("Healthy",6),
#                                 Sex=rep(NA,6),
#                                 DonorID=c("C1","C2","C1","C2","C1","C2"),
#                                 Misc=c(rep("Resting",6)))
# sample_table <- rbind(sample_table,rodriguez_samples)
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(rodriguez_cpgs),]
# rodriguez_cpgs <- rodriguez_cpgs[rownames(rodriguez_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cbind(cpg_table,rodriguez_cpgs)
# 
# #Adding Pitaksalee 2020 to the sample table.
# pitaksalee_cpgs <- preprocessDataset("Pitaksalee2020/GSE121192_series_matrix.txt")
# pitaksalee_cpgs <- pitaksalee_cpgs[,c(c(1:4),c(15:20),c(31:36))]
# pitaksalee_samples <- data.frame(ID=colnames(pitaksalee_cpgs),
#                                  Author=rep("Pitaksalee",16),
#                                  Year=rep(2020,16),
#                                  Tissue=rep("Blood",16),
#                                  CellType=c(rep("Naive CD4+",4),
#                                             rep("Memory CD4+",6),
#                                             rep("Monocyte",6)),
#                                  Age=rep(NA,16),
#                                  Condition=rep("Healthy",16),
#                                  Sex=c("F","M","F","M",
#                                        "M","F","F","M","F","M",
#                                        "M","F","F","M","F","M"),
#                                  DonorID=c("E1","E2","E3","E4",
#                                            "E5","E6","E1","E2","E3","E4",
#                                            "E5","E6", "E1","E2","E3","E4"),
#                                  Misc=c(rep("Resting",16)))
# sample_table <- rbind(sample_table,pitaksalee_samples)
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(pitaksalee_cpgs),]
# pitaksalee_cpgs <- pitaksalee_cpgs[rownames(pitaksalee_cpgs) %in% rownames(cpg_table),]
# cpg_table <- cbind(cpg_table,pitaksalee_cpgs)

#Writing tables and saving them.
#write.csv(sample_table,"ClockConstruction/validation_sample_table.csv")
#write.csv(cpg_table,"ClockConstruction/validation_cpg_table.csv")

#Reading in tables.
validation_sample_table <- read.csv("ClockConstruction/validation_sample_table.csv",row.names=1)
validation_cpg_table <- read.csv("ClockConstruction/validation_cpg_table.csv",row.names=1)
predicted_ages_newclock <- read.csv("ClockConstruction/predictions.csv",row.names=1)
validation_sample_table$new_predictions <- predicted_ages_newclock$x

#Here - testing whether or not I can re-capitulate Jonkman 2022 results.
# Looks like results are similar. 
cpg_names <- data.frame(rownames(validation_cpg_table))
validation_cpg_table <- as.data.frame(lapply(validation_cpg_table, as.numeric))
cpg_names <- cpg_names[!rowSums(validation_cpg_table > 1),]
validation_cpg_table <- validation_cpg_table[!rowSums(validation_cpg_table > 1),]
validation_cpg_table <- cbind(cpg_names,validation_cpg_table)
clocks <- checkClocks(validation_cpg_table)
predictedAges <- DNAmAge(validation_cpg_table, toBetas=FALSE)
colnames(predictedAges) <- c("ID","Horvath","Hannum","Levine","BNN",
                             "skinHorvath","PedBE","Wu","TL","BLUP","EN")
merged_data <- merge(predictedAges,validation_sample_table,by="ID")
filtered_data <- merged_data[merged_data$Misc != "Activated",]
summary <- filtered_data %>% group_by(CellType) %>% summarise(average_age = mean(new_predictions),
                                                              sd_age=sd(new_predictions))
ggplot(summary,aes(x=CellType,y=average_age)) +
  geom_bar(stat="identity",position="dodge") + geom_errorbar(
    aes(ymin=average_age-sd_age,ymax=average_age+sd_age,width=.2)) + 
  theme_classic()

