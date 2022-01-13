#load packages
library(tidyr)
library(tidyverse)

#load in MDP gtf
mdp_gtf <- read.delim("~/Desktop/Data/mdp_atg_noATC.gtf", header=FALSE)
colnames(mdp_gtf) <- c("chr","idc","idc2","start","end","idc3","sense","idc4","gene")
mdp_gtf <-separate(mdp_gtf,gene,c("gene","transcript"),";")
mdp_gtf <-separate(mdp_gtf,gene,c("trash","mdp_id"),">")
mdp_db <- mdp_gtf[,c("start","end","sense","mdp_id")]
mdp_db <- mdp_db[-c(62,425),]

#load in mitochondrial genome gtf
mitogene_gtf <- read.delim("~/Desktop/Data/mitochondria_db.gtf", header=FALSE, comment.char="#")
mitogene_db <- mitogene_gtf[,c("V3","V4","V5","V7","V9")]
mitogene_db <- mitogene_db[mitogene_db$V3 == "exon",]
colnames(mitogene_db) <- c("exon","start","end","sense","mitogene_name")
mitogene_db$mitogene_name <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1",
                               "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2",
                               "MT-TW", "MT-TA", "MT-TN", "MT-TC","MT-TY",
                               "MT-CO1", "MT-TS","MT-TD", "MT-CO2","MT-TK",
                               "MT-ATP8","MT-ATP6","MT-CO3", "MT-TG", "MT-ND3",
                               "MT-TR","MT-ND4L","MT-ND4","MT-TH","MT-TS2",
                               "MT-TL2","MT-ND5","MT-ND6","MT-TE","MT-CYB",
                               "MT-TT","MT-TP")

#set up an empty table that will hold overlap # information
encompass_table <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(encompass_table) <- c("mdp","mitogene","perc_overlap")
encompass_table$mdp <- as.character(encompass_table$mdp)
encompass_table$mitogene<- as.character(encompass_table$mitogene)
encompass_table$perc_overlap<- as.double(encompass_table$perc_overlap)

#process which MDPs are in which mitogenes and how much they overlap
for (mdp_row in 1:nrow(mdp_db)) {
  mdp_s <- mdp_db[mdp_row,"start"]
  mdp_e <- mdp_db[mdp_row,"end"]
  for (mito_row in 1:nrow(mitogene_db)) {
    mito_s <- mitogene_db[mito_row,"start"]
    mito_e <- mitogene_db[mito_row,"end"]
    if (
      (mdp_s < mito_s & mdp_s < mito_e & mdp_e < mito_s & mdp_e < mito_e) | 
      (mdp_s > mito_s & mdp_s > mito_e & mdp_e > mito_s & mdp_e > mito_e)) {
      next;
    }
    if (mdp_s > mito_s & mdp_s < mito_e & mdp_e < mito_s & mdp_e < mito_e) {
      overlap <- (mdp_e - mito_s) / (mdp_e - mdp_s)
      encompass_table %>% add_row(mdp=mdp_db[mdp_row,"mdp_id"],
                                  mitogene = mitogene_db[mito_row,"mitogene_name"],
                                  perc_overlap = overlap)
      break;
    }
    if (mdp_s > mito_s & mdp_s < mito_e & mdp_e > mito_s & mdp_e > mito_e) {
      overlap <- (mito_e - mdp_s) / (mdp_e - mdp_s)
      encompass_table %>% add_row(mdp=mdp_db[mdp_row,"mdp_id"],
                                  mitogene = mitogene_db[mito_row,"mitogene_name"],
                                  perc_overlap = overlap)
      break;
    }
    if (mdp_s > mito_s & mdp_s < mito_e & mdp_e > mito_s & mdp_e < mito_e) {
      overlap <- 1
      encompass_table <- encompass_table %>% 
        add_row(mdp = mdp_db[mdp_row,"mdp_id"],
                mitogene = mitogene_db[mito_row,"mitogene_name"],
                perc_overlap = overlap)
      break;
    }
  } 
}


