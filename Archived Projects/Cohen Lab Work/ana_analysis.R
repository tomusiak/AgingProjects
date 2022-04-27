flowpbmcs <- read.csv("~/Desktop/Data/ana_data.csv")
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")
library(ggplot2)
library(RColorBrewer)
library(dplyr)

data_summary_activatedcd8s <- getSummary(flowpbmcs,"CD8_Activated","Condition")
ggplot(data=data_summary_activatedcd8s,aes(x=Condition,y=CD8_Activated)) + 
  geom_bar(stat="identity",position="dodge",fill="maroon",color="black") +
  geom_errorbar(aes(ymin=CD8_Activated-se, ymax=CD8_Activated+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="% CD69+ (across all CD8s+)",title="MDP Effect on Activation of CD8+ Cells") +
  theme_classic() + scale_fill_brewer(palette="Paired")

data_summary_activatedcd4s <- getSummary(flowpbmcs,"CD4_Activated","Condition")
ggplot(data=data_summary_activatedcd4s,aes(x=Condition,y=CD4_Activated)) + 
  geom_bar(stat="identity",position="dodge",fill="maroon",color="black") +
  geom_errorbar(aes(ymin=CD4_Activated-se, ymax=CD4_Activated+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="% CD69+ (across all CD8s+)",title="MDP Effect on Activation of CD8+ Cells") +
  theme_classic() + scale_fill_brewer(palette="Paired")

flowpbmcs$Cell_Count <- flowpbmcs$Living.Cells * (10000/flowpbmcs$Beads)

data_summary_count <- getSummary(flowpbmcs,"Cell_Count","MDP")
data_summary_count$MDP <-factor(data_summary_count$MDP,levels=c("No MDP","151C","18A","159C",
                                                                "6A","NOSH","140B","144B"))
ggplot(data=data_summary_count,aes(x=MDP,y=Cell_Count)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=Cell_Count-se, ymax=Cell_Count+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Total Alive Cells",title="MDP Effects on Cell Counts") +
  theme_classic() + scale_fill_brewer(palette="Paired")

data_summary_naivecd8 <- getSummary(flowpbmcs,"CD8_Naive","MDP")
data_summary_naivecd8$MDP <-factor(data_summary_naivecd8$MDP,levels=c("151C","18A","159C",
                                                                      "6A","NOSH","140B","144B","No MDP"))
ggplot(data=data_summary_naivecd8,aes(x=MDP,y=CD8_Naive)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=CD8_Naive-se, ymax=CD8_Naive+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Total Alive Cells",title="%Naive (Relative to CD8+)") +
  theme_classic() + scale_fill_brewer(palette="Paired")
data_summary_naivecd8$celltype <- "Naive"
colnames(data_summary_naivecd8) <- c("MDP","N","Percent","SD","SE","CI","celltype")

data_summary_cmcd8 <- getSummary(flowpbmcs,"CD8_CM","MDP")
data_summary_cmcd8$MDP <-factor(data_summary_cmcd8$MDP,levels=c("151C","18A","159C",
                                                                "6A","NOSH","140B","144B","No MDP"))
ggplot(data=data_summary_cmcd8,aes(x=MDP,y=CD8_CM)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=CD8_CM-se, ymax=CD8_CM+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Total Alive Cells",title="%CM (Relative to CD8+)") +
  theme_classic() + scale_fill_brewer(palette="Paired")
data_summary_cmcd8$celltype <- "CM"
colnames(data_summary_cmcd8) <- c("MDP","N","Percent","SD","SE","CI","celltype")

data_summary_emcd8 <- getSummary(flowpbmcs,"CD8_EM","MDP")
data_summary_emcd8$MDP <-factor(data_summary_emcd8$MDP,levels=c("151C","18A","159C",
                                                                "6A","NOSH","140B","144B","No MDP"))
ggplot(data=data_summary_emcd8,aes(x=MDP,y=CD8_EM)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=CD8_EM-se, ymax=CD8_EM+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Total Alive Cells",title="%CM (Relative to CD8+)") +
  theme_classic() + scale_fill_brewer(palette="Paired")
data_summary_emcd8$celltype <- "EM"
colnames(data_summary_emcd8) <- c("MDP","N","Percent","SD","SE","CI","celltype")

data_summary_temracd8 <- getSummary(flowpbmcs,"CD8_TEMRA","MDP")
data_summary_temracd8$MDP <-factor(data_summary_temracd8$MDP,levels=c("151C","18A","159C",
                                                                      "6A","NOSH","140B","144B","No MDP"))
ggplot(data=data_summary_temracd8,aes(x=MDP,y=CD8_TEMRA)) + 
  geom_bar(stat="identity",position="dodge",fill="lightblue",color="black") +
  geom_errorbar(aes(ymin=CD8_TEMRA-se, ymax=CD8_TEMRA+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="Total Alive Cells",title="%CM (Relative to CD8+)") +
  theme_classic() + scale_fill_brewer(palette="Paired")
data_summary_temracd8$celltype <- "TEMRA"
colnames(data_summary_temracd8) <- c("MDP","N","Percent","SD","SE","CI","celltype")

all_subsets <- rbind(data_summary_naivecd8,data_summary_cmcd8,data_summary_emcd8,data_summary_temracd8)
ggplot(data=all_subsets,aes(x=MDP,y=Percent,fill=celltype)) + 
  geom_bar(stat="identity",position="dodge") +
  geom_errorbar(aes(ymin=Percent-SE, ymax=Percent+SE), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (10uM)",y="% (of total CD8+)",title="CD8+ Cell Type Composition") +
  theme_classic() + scale_fill_brewer(palette="Paired")

