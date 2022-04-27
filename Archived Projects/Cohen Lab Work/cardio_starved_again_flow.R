flowcardiomyocytes <- read.csv("~/Desktop/Data/signal.csv")
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")
library(ggplot2)
library(RColorBrewer)
library(dplyr)

starved <- flowcardiomyocytes

starved_alive <- getSummary(starved,"Alive_Percent","MDP")
starved_alive$MDP <- factor(starved_alive$MDP,levels=c("No Doxo","Doxo Only","63D WT","63D SNP",
                                                       "108C WT","108C SNP","HNG","NOSH"))
ggplot(data=starved_alive,aes(x=MDP,y=Alive_Percent)) + 
  geom_bar(stat="identity",position="dodge",fill="pink",color="black") +
  geom_errorbar(aes(ymin=Alive_Percent-se, ymax=Alive_Percent+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="% Alive",title="MDP Effects on Doxo-Induced Cell Viability") +
  theme_classic() + scale_fill_brewer(palette="Paired")

starved_count <- getSummary(starved,"Cells","MDP")
starved_count$MDP <- factor(starved_count$MDP,levels=c("No Doxo","Doxo Only","63D WT","63D SNP",
                                                       "108C WT","108C SNP","HNG","NOSH"))
ggplot(data=starved_count,aes(x=MDP,y=Cells)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=Cells-se, ymax=Cells+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Cell-like Objects Per 100uL",title="MDP Protection Against Doxo-Induced Changes in Cell-like Objects") +
  theme_classic() + scale_fill_brewer(palette="Paired")

starved_roshigh <- getSummary(starved,"ROS_Percent","MDP")
starved_roshigh$MDP <- factor(starved_roshigh$MDP,levels=c("No Doxo","Doxo Only","63D WT","63D SNP",
                                                           "108C WT","108C SNP","HNG","NOSH"))
ggplot(data=starved_roshigh,aes(x=MDP,y=ROS_Percent)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=ROS_Percent-se, ymax=ROS_Percent+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="ROS+ % Cells",title="MDP Protection Against Doxo-Induced ROS Production") +
  theme_classic() + scale_fill_brewer(palette="Paired")

starved_roshigh_signal <- getSummary(starved,"ROS_Signal","MDP")
starved_roshigh_signal$MDP <- factor(starved_roshigh_signal$MDP,levels=c("No Doxo","Doxo Only","63D WT","63D SNP",
                                                           "108C WT","108C SNP","HNG","NOSH"))
ggplot(data=starved_roshigh_signal,aes(x=MDP,y=ROS_Signal)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=ROS_Signal-se, ymax=ROS_Signal+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="ROS+ % Cells",title="MDP Protection Against Doxo-Induced ROS Production") +
  theme_classic() + scale_fill_brewer(palette="Paired")

alive_signal <- getSummary(starved,"Calcein_Signal","MDP")
alive_signal$MDP <- factor(alive_signal$MDP,levels=c("No Doxo","Doxo Only","63D WT","63D SNP",
                                                                         "108C WT","108C SNP","HNG","NOSH"))
ggplot(data=alive_signal,aes(x=MDP,y=Calcein_Signal)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=Calcein_Signal-se, ymax=Calcein_Signal+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Alive Signal",title="MDP Effects on Doxo-Induced Cell Viability (Mean Alive Signal)") +
  theme_classic() + scale_fill_brewer(palette="Paired")
