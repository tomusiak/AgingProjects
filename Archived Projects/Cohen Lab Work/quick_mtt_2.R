mtt_data <-read.csv("~/Desktop/AgingProjects/Cohen/3-31mtt.csv")
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")
library(ggplot2)
library(RColorBrewer)
library(dplyr)

data_summary <- getSummary(mtt_data,"Values","MDP")
data_summary$MDP <-factor(data_summary$MDP,levels=c("No Doxo",
                  "63D WT", 
                  "63D SNP", 
                  "108C WT", 
                  "108C SNP",
                  "HNG", "NOSH",
                  "Doxo Only"))

ggplot(data=data_summary,aes(x=MDP,y=Values)) + 
  geom_bar(stat="identity",position="dodge",fill="pink",color="black") +
  geom_errorbar(aes(ymin=Values-se, ymax=Values+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide (100uM)",y="Absorbance (Ab550-Ab680)",title="MDP Protection Against Doxo Toxicity") +
  theme_classic() + scale_fill_brewer(palette="Paired") + coord_flip()

