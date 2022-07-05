#Grabs some useful scripts.
source("AgingProjects/Useful Scripts/generally_useful.R")

#Sets location of data 
setwd("Data/") #Sets directory.
#Libraries to import.
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(plyr)
library(reshape2)
library(WGCNA)
library(dplyr)
library("RColorBrewer")
library(readxl)
library(ggplot2)
library(GOfuncR)
library(readr)
library(umap)
library(stringr)
library(tidyr)
library(gplots)
library(minfi)

#Read in clock data and metadata. Merge them and process them.
clock_data <- read.csv("subset_data/clock_data.csv")
subset_metadata <- read.csv("subset_data/subset_metadata.csv")
all_data <- merge(clock_data,subset_metadata)
all_data$type <- as.character(all_data$type)
all_data$type <- factor(all_data$type, levels=c("naive", "central_memory", "effector_memory","temra"))

beta_values <- read.csv("subset_data/beta_values.csv", row.names=1)
m_values <- read.csv("subset_data/m_values.csv", row.names=1)

#Filter out unwanted data from metadata, mapping, and beta values.
keep <- all_data$SampleID[all_data$SampleID != "D3" & all_data$sabgal_sample == FALSE]
all_data <- all_data[all_data$SampleID %in% keep,]
beta_values <- beta_values[,colnames(beta_values) %in% keep]
beta_values <- beta_values[order(all_data$type)]
m_values <- m_values[,colnames(m_values) %in% keep]
m_values <- m_values[order(all_data$type)]
all_data <- all_data[order(all_data$type),]

#Calculate epigenetic clock acceleration
all_data$diff <- all_data$DNAmAge-all_data$age
summary <- getSummary(all_data,"diff", "type")
summary$type <- c("Naive","Central Memory","Effector Memory","TEMRA")
summary$type <- factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

#subset ages
naive_data <- all_data[all_data$type == "naive",]
cm_data <- all_data[all_data$type == "central_memory",]
em_data <- all_data[all_data$type == "effector_memory",]
temra_data <- all_data[all_data$type == "temra",]

#t tests
t.test(em_data$IEAA,naive_data$IEAA,paired=TRUE)
t_test_data <- data.frame(all_data[all_data$donor != "R45553",])
t_test_data <- t_test_data[c("IEAA","type","donor")]
t_test_data$type <- as.factor(t_test_data$type)
t_test_data$donor <- as.factor(t_test_data$donor)
test <- t_test_data %>% pairwise_t_test(IEAA ~ type, paired=TRUE,correction="bonferroni")  %>%
  select(-df, -statistic, -p)
res.aov <- t_test_data %>% anova_test(IEAA ~ type) 
res.aov

#QC check for DNA quantity
ggplot(data=all_data, aes(x=as.factor(final_dna), y=meanAbsDifferenceSampleVSgoldstandard, group=1)) +
  geom_point() +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  labs(x="Input DNA",y="Difference between Sample and Gold Standard",
       title="QC - Differences between Sample and Gold Standard") 

ggplot(data=summary, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-20,20) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age") 

summary_IEAA <- getSummary(all_data,"IEAA", "type")

ggplot(data=summary_IEAA, aes(x=type, y=IEAA, group=1)) +
  geom_point() +
  theme_classic() +
  ylim(-5,5) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=IEAA-se, ymax=IEAA+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="IEAA", title="CD8+ T Cell Subset Differences With IEAA") 

summary_IEAA_donor <- getSummary(all_data,"IEAA", "donor")

ggplot(data=all_data, aes(x=donor, y=IEAA, group=1,color=type)) +
  geom_point() +
  theme_classic() +
  ylim(-10,10) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  labs(x="Donor",y="IEAA", title="CD8+ IEAA Differences Between Donors") 

young_subset <- all_data[all_data$age <= 50,]
old_subset <- all_data[all_data$age > 50,]

summary_young <- getSummary(young_subset,"diff", "type")
summary_young$type <-  c("Naive","Central Memory","Effector Memory","TEMRA")
summary_young$type <-factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=summary_young, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-23,23) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age - < 50 years old") 

summary_old <- getSummary(old_subset,"diff", "type")
summary_old$type <- c("Naive","Central Memory","Effector Memory","TEMRA")
summary_old$type <- factor(summary$type,levels=c("Naive","Central Memory","Effector Memory","TEMRA"))

ggplot(data=summary_old, aes(x=type, y=diff, group=1)) +
  geom_line()+ 
  geom_point() +
  theme_classic() +
  ylim(-23,23) +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_errorbar(aes(ymin=diff-se, ymax=diff+se), width=.1) +
  labs(x="CD8+ T Cell Subset",y="Predicted Age - Age", title="CD8+ T Cell Subset Differences Between
       Clock Age and Chronological Age - >50 years old") 

ggplot(data=all_data,aes(x=age, y=DNAmAge, color=type)) +
  geom_point() +
  theme_classic() +
  xlim(10,80) + ylim(10,80) +
  theme(text = element_text(size = 15)) +
  geom_abline(linetype="dotted") +
  geom_hline(yintercept = 0,linetype="dotted") +
  labs(x="Donor Age",y="Predicted Age", title="Age vs. Predicted Age Per Donor") 

#Looking at clock CpGs
#data(HorvathLongCGlist)
#clock_list <- HorvathLongCGlist
#matching_cpgs <- match(clock_list$MR_var,mapping$cpg)
#gene_names <- mapping$row_name[matching_cpgs]
#clock_changes <- na.omit(diff_exp[gene_names,1:6])
#clock_changes <- clock_changes[order(clock_changes$adj.P.Val),]
#clock_changes$color=factor(case_when(clock_changes$adj.P.Val < .05 & abs(clock_changes$logFC) >= .3 ~ "purple",
#                                     (clock_changes$adj.P.Val < .05 & abs(clock_changes$logFC) < .3) ~ "red",
#                                     (clock_changes$adj.P.Val > .05 & abs(clock_changes$logFC) >= .3) ~ "blue",
#                                     (clock_changes$adj.P.Val > .05 & abs(clock_changes$logFC) < .3) ~ "gray"))
#clock_changes$delabel <- NA
#clock_changes$delabel[clock_changes$color=="purple"] <- rownames(clock_changes)[clock_changes$color=="purple"]
#ggplot(data=clock_changes, aes(x=logFC, y=-log10(adj.P.Val), color=color,label=delabel)) + 
#  geom_point() +
#  xlim(-1,1) + ylim(0,12) +
#  theme_classic(base_size=15)  + 
#  geom_vline(xintercept=.25,linetype="dotted") +
#  geom_vline(xintercept=-.25,linetype="dotted") +
#  geom_hline(yintercept = 1.2,linetype="dotted") +
#  geom_text(nudge_y=.2) +
#  scale_colour_identity() +
#  labs(title="Volcano Plot of Clock CpG Sites Changing between Naive & TEMRA Samples")
#head(clock_changes)