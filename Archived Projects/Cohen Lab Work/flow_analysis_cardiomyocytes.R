flowcardiomyocytes <- read.csv("~/Desktop/Data/flowcardiomyocytes.csv")
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
colnames(flowcardiomyocytes) <- c("MDP","Counting_Beads","Alive_Percent","KI67_Percent","P21_Percent",
                                  "P21_Mean_Signal","Alive_Mean_Signal","Alive_Cells","Replicate")
data_summary_alive <- getSummary(flowcardiomyocytes,"Alive_Percent","MDP")
data_summary_alive$MDP <- factor(data_summary_alive$MDP,levels=c("No Doxo","63D WT","63D SNP",
                                                                    "108C WT","108C SNP","HNG","NOSH","Doxo Only"))
ggplot(data=data_summary_alive,aes(x=MDP,y=Alive_Percent)) + 
  geom_bar(stat="identity",position="dodge",fill="pink",color="black") +
  geom_errorbar(aes(ymin=Alive_Percent-se, ymax=Alive_Percent+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="% Not Cell Permeable",title="MDP Protection Against Doxo-Induced Cell Permeability") +
  theme_classic() + scale_fill_brewer(palette="Paired")

data_summary_prolif <- getSummary(flowcardiomyocytes,"KI67_Percent","MDP")
data_summary_prolif$MDP <- factor(data_summary_alive$MDP,levels=c("No Doxo","63D WT","63D SNP",
                                                                 "108C WT","108C SNP","HNG","NOSH","Doxo Only"))
ggplot(data=data_summary_prolif,aes(x=MDP,y=KI67_Percent)) + 
  geom_bar(stat="identity",position="dodge",fill="maroon",color="black") +
  geom_errorbar(aes(ymin=KI67_Percent-se, ymax=KI67_Percent+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="% Ki-67+",title="MDP Protection Against Doxo-Induced Changes in Ki-67%") +
  theme_classic() + scale_fill_brewer(palette="Paired")

data_summary_p21 <- getSummary(flowcardiomyocytes,"P21_Percent","MDP")
data_summary_p21$MDP <- factor(data_summary_p21$MDP,levels=c("No Doxo","63D WT","63D SNP",
                                                                  "108C WT","108C SNP","HNG","NOSH","Doxo Only"))
ggplot(data=data_summary_p21,aes(x=MDP,y=P21_Percent)) + 
  geom_bar(stat="identity",position="dodge",fill="violet",color="black") +
  geom_errorbar(aes(ymin=P21_Percent-se, ymax=P21_Percent+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="% P-21+",title="MDP Protection Against Doxo-Induced Changes in Cell Cycle Arrest") +
  theme_classic() + scale_fill_brewer(palette="Paired")

data_summary_count <- getSummary(flowcardiomyocytes,"Alive_Cells","MDP")
data_summary_count$MDP <- factor(data_summary_count$MDP,levels=c("No Doxo","63D WT","63D SNP",
                                                               "108C WT","108C SNP","HNG","NOSH","Doxo Only"))
ggplot(data=data_summary_count,aes(x=MDP,y=Alive_Cells)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=Alive_Cells-se, ymax=Alive_Cells+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Alive Cells Per 10uL",title="MDP Protection Against Doxo-Induced Changes in Cell Counts") +
  theme_classic() + scale_fill_brewer(palette="Paired")
