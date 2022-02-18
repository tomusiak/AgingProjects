mtt_data <- read.csv("~/Desktop/Data/2_17_22_mtt.csv")
setwd("/home/atom/Desktop/AgingProjects/Useful Scripts/")
source("generally_useful.R")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
mtt_data$Combo <- paste(mtt_data$concentration_uM,mtt_data$mdp)
data_summary <- getSummary(mtt_data,"value","Combo")

mtt_data <- mtt_data %>% group_by(mdp,concentration_uM) %>% summarize(mean_value = mean(value),
                                                                          sd = sd(value),
                                                                          n = n(),
                                                                          se = sd / sqrt(n))
mtt_data$mdp <-factor(mtt_data$mdp,levels=c("108C SNP",
                  "108C WT", 
                  "63D SNP", 
                  "63D WT", 
                  "MOTS_C",
                  "Doxo-", "Doxo+"))

mtt_data$mdp <- as.factor(mtt_data$mdp)

ggplot(data=mtt_data,aes(x=as.factor(concentration_uM),y=mean_value, fill=mdp)) + 
  geom_bar(stat="identity",position="dodge") +
  geom_errorbar(aes(ymin=mean_value-se, ymax=mean_value+se), width=.2,position=position_dodge(width=.9)) +
  labs(x="Peptide Concentration (uM)",y="Absorbance (Ab550-Ab680)",title="MDP Protection Against Doxo Toxicity") +
  theme_classic() + scale_fill_brewer(palette="Paired")

