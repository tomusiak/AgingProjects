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

# #Quick helper function for pre-processing datasets.
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
# validation_sample_table <- rbind(sample_table,pitaksalee_samples)
# cpg_table <- cpg_table[rownames(cpg_table) %in% rownames(pitaksalee_cpgs),]
# pitaksalee_cpgs <- pitaksalee_cpgs[rownames(pitaksalee_cpgs) %in% rownames(cpg_table),]
# validation_cpg_table <- cbind(cpg_table,pitaksalee_cpgs)
# 
# #Writing tables and saving them.
# write.csv(validation_sample_table,"ClockConstruction/validation_sample_table.csv")
# write.csv(validation_cpg_table,"ClockConstruction/validation_cpg_table.csv")

#Reading in tables.
validation_sample_table <- read.csv("ClockConstruction/validation_sample_table.csv",row.names=1)
validation_cpg_table <- read.csv("ClockConstruction/validation_cpg_table.csv",row.names=1)

# #Here - testing whether or not I can re-capitulate Jonkman 2022 results.
# # Looks like results are similar. 
# cpg_names <- data.frame(rownames(validation_cpg_table))
# validation_cpg_table <- as.data.frame(lapply(validation_cpg_table, as.numeric))
# cpg_names <- cpg_names[!rowSums(validation_cpg_table > 1),]
# validation_cpg_table <- validation_cpg_table[!rowSums(validation_cpg_table > 1),]
# validation_cpg_table <- cbind(cpg_names,validation_cpg_table)
# clocks <- checkClocks(validation_cpg_table)
# predicted_ages <- DNAmAge(validation_cpg_table, toBetas=FALSE)
# colnames(predicted_ages) <- c("ID","Horvath","Hannum","Levine","BNN",
#                              "skinHorvath","PedBE","Wu","TL","BLUP","EN")
# write.csv(predicted_ages,"ClockConstruction/predicted_ages_oldclock.csv")

predicted_ages_oldclock <- read.csv("ClockConstruction/predicted_ages_oldclock.csv",row.names=1)
predicted_ages_newclock <- read.csv("ClockConstruction/predictions.csv",row.names=1)
validation_sample_table$new_predictions <- predicted_ages_newclock$x
predicted_ages_newclock$ID <- rownames(predicted_ages_newclock)
merged_data <- merge(predicted_ages_newclock,validation_sample_table,by="ID")
merged_data <- merge(merged_data,predicted_ages_oldclock,by="ID")
filtered_data <- merged_data[merged_data$Misc != "Activated",]

#Going to investigate CD4s and CD8s and compare the "differences" between CD8+ effector -> naive
# and CD4 memory -> naive. If the new clock has differences substantially closer to 0, that is
# an indicator of success.
cd8_summary_type <- filtered_data[filtered_data$CellType=="Effector CD8+" | 
                                filtered_data$CellType=="Naive CD8+",] %>% 
  group_by(DonorID,CellType) %>% summarise(average_age = mean(Horvath))
cd8_differences <- cd8_summary_type %>% group_by(DonorID)
ggplot(cd8_differences,aes(x=CellType,y=average_age)) + 
  geom_line(aes(group=DonorID)) + theme_classic()

cd8_new_summary_type <- filtered_data[filtered_data$CellType=="Effector CD8+" | 
                                filtered_data$CellType=="Naive CD8+",] %>% 
  group_by(DonorID,CellType) %>% summarise(average_age = mean(new_predictions))
cd8_new_differences <- cd8_new_summary_type %>% group_by(DonorID)
ggplot(cd8_new_summary_type,aes(x=CellType,y=average_age)) + 
  geom_line(aes(group=DonorID)) + theme_classic()

cd8_old_changes<- na.omit(ave(cd8_differences$average_age, 
                              factor(cd8_differences$DonorID), FUN=function(x) c(NA,diff(x))))
cd8_new_changes<- na.omit(ave(cd8_new_differences$average_age, 
                              factor(cd8_new_differences$DonorID), FUN=function(x) c(NA,diff(x))))
cd8_old_mean <- mean(cd8_old_changes)
cd8_old_se <- sd(cd8_old_changes)/sqrt(6)
cd8_new_mean <- mean(cd8_new_changes)
cd8_new_se <- sd(cd8_new_changes)/sqrt(6)

cd8_quantifying_differences <- data.frame("Mean_Difference" = 1,
                                      "Mean_SE" = 1)

cd8_quantifying_differences <- rbind(cd8_quantifying_differences,
                                     data.frame(Mean_Difference=cd8_old_mean,
                                                Mean_SE=cd8_old_se))
cd8_quantifying_differences <- rbind(cd8_quantifying_differences,
                                     data.frame(Mean_Difference=cd8_new_mean,
                                                Mean_SE=cd8_new_se))
cd8_quantifying_differences <- cd8_quantifying_differences[-1,]
cd8_quantifying_differences$Clock <- c("Horvath_Clock", "New_Clock")

#Let's do a T test with the old and new clock on the validation samples - CD8s.
old_naive_cd8 <- cd8_summary_type[cd8_summary_type$CellType=="Naive CD8+",]
old_effector_cd8 <- cd8_summary_type[cd8_summary_type$CellType=="Effector CD8+",]
new_naive_cd8 <- cd8_new_summary_type[cd8_new_summary_type$CellType=="Naive CD8+",]
new_effector_cd8 <- cd8_new_summary_type[cd8_new_summary_type$CellType=="Effector CD8+",]
t.test(old_naive_cd8$average_age,old_effector_cd8$average_age,paired=TRUE)
t.test(new_naive_cd8$average_age,new_effector_cd8$average_age,paired=TRUE)

#CD8+ differentiation processed above. Will process CD4+ differentiation below.

cd4_summary_type <- filtered_data[filtered_data$CellType=="Naive CD4+" | 
                                    filtered_data$CellType=="Memory CD4+",] %>% 
  group_by(DonorID,CellType) %>% summarise(average_age = mean(Horvath))
cd4_summary_type <- cd4_summary_type[1:14,]
cd4_differences <- cd4_summary_type %>% group_by(DonorID)
ggplot(cd4_differences,aes(x=CellType,y=average_age)) + 
  geom_line(aes(group=DonorID)) + theme_classic()

cd4_new_summary_type <- filtered_data[filtered_data$CellType=="Naive CD4+" | 
                                    filtered_data$CellType=="Memory CD4+",] %>% 
  group_by(DonorID,CellType) %>% summarise(average_age = mean(new_predictions))
cd4_new_summary_type <- cd4_new_summary_type[1:14,]
cd4_new_differences <- cd4_new_summary_type %>% group_by(DonorID)
ggplot(cd4_new_differences,aes(x=CellType,y=average_age)) + 
  geom_line(aes(group=DonorID)) + theme_classic()

cd4_old_changes<- na.omit(ave(cd4_differences$average_age, 
                              factor(cd4_differences$DonorID), FUN=function(x) c(NA,diff(x))))
cd4_new_changes<- na.omit(ave(cd4_new_differences$average_age, 
                              factor(cd4_new_differences$DonorID), FUN=function(x) c(NA,diff(x))))
cd4_old_mean <- mean(cd4_old_changes)
cd4_old_se <- sd(cd4_old_changes)/sqrt(7)
cd4_new_mean <- mean(cd4_new_changes)
cd4_new_se <- sd(cd4_new_changes)/sqrt(7)

cd4_quantifying_differences <- data.frame("Mean_Difference" = 1,
                                      "Mean_SE" = 1)

cd4_quantifying_differences <- rbind(cd4_quantifying_differences,
                                     data.frame(Mean_Difference=cd4_old_mean,
                                                Mean_SE=cd4_old_se))
cd4_quantifying_differences <- rbind(cd4_quantifying_differences,
                                     data.frame(Mean_Difference=cd4_new_mean,
                                                Mean_SE=cd4_new_se))
cd4_quantifying_differences <- cd4_quantifying_differences[-1,]
cd4_quantifying_differences$Clock <- c("Horvath_Clock", "New_Clock")
       
cd8_quantifying_differences$Cell <- "CD8 EM - CD8 Naive"
cd4_quantifying_differences$Cell <- "CD4 Memory - CD4 Naive"

quantifying_differences <- rbind(cd8_quantifying_differences,
                                 cd4_quantifying_differences)

#Let's do a T test with the old and new clock on the validation samples - CD4s.
old_naive_cd4 <- cd4_summary_type[cd4_summary_type$CellType=="Naive CD4+",]
old_memory_cd4 <- cd4_summary_type[cd4_summary_type$CellType=="Memory CD4+",]
new_naive_cd4 <- cd4_new_summary_type[cd4_new_summary_type$CellType=="Naive CD4+",]
new_memory_cd4 <- cd4_new_summary_type[cd4_new_summary_type$CellType=="Memory CD4+",]
t.test(old_naive_cd4$average_age,old_memory_cd4$average_age,paired=TRUE)
t.test(new_naive_cd4$average_age,new_memory_cd4$average_age,paired=TRUE)

# Now plotting the differences in predicted age per differentiation state per clock.

ggplot(quantifying_differences, aes(x=Clock,y=Mean_Difference,
                                        ymin=Mean_Difference-Mean_SE,ymax=Mean_Difference+Mean_SE)) +
  geom_point() + theme_classic() + geom_errorbar() + ylim(-25,25) + facet_wrap(~Cell) +
  geom_hline(yintercept=0,linetype="dashed")

