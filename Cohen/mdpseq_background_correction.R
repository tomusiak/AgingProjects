#load packages
library(tidyr)
library(tidyverse)

read_length <- 90

generateEncompassTable <- function(mdp_gtf, mitochondrial_gtf) {
  #Pre-processing MDP gtf.
  colnames(mdp_gtf) <-
    c("chr",
      "idc",
      "idc2",
      "start",
      "end",
      "idc3",
      "sense",
      "idc4",
      "gene")
  mdp_gtf <- separate(mdp_gtf, gene, c("gene", "transcript"), ";")
  mdp_gtf <- separate(mdp_gtf, gene, c("trash", "mdp_id"), ">")
  mdp_db <- mdp_gtf[, c("start", "end", "sense", "mdp_id")]
  mdp_db <- mdp_db[-c(62, 425), ]
  
  #Pre-processing mitochondrial GTF.
  mitogene_db <- mitogene_gtf[, c("V3", "V4", "V5", "V7", "V9")]
  mitogene_db <- mitogene_db[mitogene_db$V3 == "exon", ]
  colnames(mitogene_db) <-
    c("exon", "start", "end", "sense", "mitogene_name")
  mitogene_db$mitogene_name <-
    c(
      "MT-TF",
      "MT-RNR1",
      "MT-TV",
      "MT-RNR2",
      "MT-TL1",
      "MT-ND1",
      "MT-TI",
      "MT-TQ",
      "MT-TM",
      "MT-ND2",
      "MT-TW",
      "MT-TA",
      "MT-TN",
      "MT-TC",
      "MT-TY",
      "MT-CO1",
      "MT-TS",
      "MT-TD",
      "MT-CO2",
      "MT-TK",
      "MT-ATP8",
      "MT-ATP6",
      "MT-CO3",
      "MT-TG",
      "MT-ND3",
      "MT-TR",
      "MT-ND4L",
      "MT-ND4",
      "MT-TH",
      "MT-TS2",
      "MT-TL2",
      "MT-ND5",
      "MT-ND6",
      "MT-TE",
      "MT-CYB",
      "MT-TT",
      "MT-TP"
    )
  
  #Setting up an empty table that will hold overlap # information
  encompass_table <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(encompass_table) <- c("mdp", "mitogene", "perc_overlap","proportion")
  encompass_table$mdp <- as.character(encompass_table$mdp)
  encompass_table$mitogene <- as.character(encompass_table$mitogene)
  encompass_table$perc_overlap <-
    as.double(encompass_table$perc_overlap)
  encompass_table$proportion <-
    as.double(encompass_table$proportion)
  
  #Processing which MDPs are in which mitogenes and how much they overlap
  for (mdp_row in 1:nrow(mdp_db)) {
    mdp_s <- mdp_db[mdp_row, "start"]
    mdp_e <- mdp_db[mdp_row, "end"]
    mdp_size <- mdp_e-mdp_s
    for (mito_row in 1:nrow(mitogene_db)) {
      mito_s <- mitogene_db[mito_row, "start"]
      mito_e <- mitogene_db[mito_row, "end"]
      mito_size <- mito_e-mito_s
      if ((mdp_s < mito_s &
           mdp_s < mito_e & mdp_e < mito_s & mdp_e < mito_e) |
          (mdp_s > mito_s &
           mdp_s > mito_e & mdp_e > mito_s & mdp_e > mito_e)) {
        next;
      }
      if (mdp_s < mito_s &
          mdp_s < mito_e & mdp_e > mito_s & mdp_e < mito_e) {
        overlap <- (mdp_e - mito_s) / (mdp_e - mdp_s)
        encompass_table <- encompass_table %>% add_row(mdp = mdp_db[mdp_row, "mdp_id"],
                                    mitogene = mitogene_db[mito_row, "mitogene_name"],
                                    perc_overlap = overlap, proportion=overlap*((mdp_size+read_length)/mito_size))
      } 
      if (mdp_s > mito_s &
          mdp_s < mito_e & mdp_e > mito_s & mdp_e > mito_e) {
        overlap <- (mito_e - mdp_s) / (mdp_e - mdp_s)
        encompass_table <- encompass_table %>% add_row(mdp = mdp_db[mdp_row, "mdp_id"],
                                    mitogene = mitogene_db[mito_row, "mitogene_name"],
                                    perc_overlap = overlap,proportion=overlap*((mdp_size+read_length)/mito_size))
      } 
      if (mdp_s > mito_s &
          mdp_s < mito_e & mdp_e > mito_s & mdp_e < mito_e) {
        overlap <- 1
        encompass_table <- encompass_table %>%
          add_row(mdp = mdp_db[mdp_row, "mdp_id"],
                  mitogene = mitogene_db[mito_row, "mitogene_name"],
                  perc_overlap = overlap,
                  proportion = (mdp_size+read_length)/mito_size)
      } 
    }
  }
  
  for (mdp_row in 1:nrow(encompass_table)) {
    outer_mdp <- encompass_table$mdp[mdp_row]
    outer_overlap <- encompass_table$perc_overlap[mdp_row]
    for (mdp_row_2 in mdp_row+1:nrow(encompass_table)) {
      inner_mdp <- encompass_table$mdp[mdp_row_2]
      inner_overlap <- encompass_table$perc_overlap[mdp_row_2]
      if (!is.na(inner_mdp)) {
        if (outer_mdp == inner_mdp) {
          if (outer_overlap > inner_overlap) {
            encompass_table <- encompass_table[-mdp_row_2,]
          }
          if (outer_overlap == inner_overlap) {
            encompass_table <- encompass_table[-mdp_row_2,]
          }
          if (outer_overlap < inner_overlap) {
            encompass_table <- encompass_table[-mdp_row,]
          }
        }
      }
    }
  }
  
  #Spitting out the data.
  return (encompass_table)
}

determineBackgroundSignal <- function(mito_counts) {
  average <- rowMeans(mito_counts[1:ncol(mito_counts)], na.rm=TRUE)
  mito_counts$average <- average
  for (mito_column in 1:(ncol(mito_counts)-1)) {
      mito_counts[,mito_column] <- (mito_counts[,mito_column] - mito_counts[,ncol(mito_counts)]) /
        (mito_counts[,ncol(mito_counts)])
  }
  mito_counts <- mito_counts[,-ncol(mito_counts)]
  return (mito_counts)
}

performBackgroundCorrection <- function(background_table,mdp_counts,encompass_table) {
  mdp_names <- rownames(mdp_counts)
  for (mdp_row in 1:nrow(mdp_counts)) {
    mdp <- mdp_names[mdp_row]
    corres_mitogene <- encompass_table$mitogene[encompass_table$mdp==mdp]
    corres_overlap <- encompass_table$perc_overlap[encompass_table$mdp==mdp]
    corres_proportion <- 1-encompass_table$proportion[encompass_table$mdp==mdp]
    if (length(corres_proportion) != 0) {
      if (corres_proportion < 0) {
        corres_proportion <- 0
      } 
    }
    if(length(corres_mitogene) != 0) {
      for (mdp_col in 2:length(mdp_counts[mdp_row,])) {
        initial_mdp_count <- mdp_counts[mdp_row, 1]
        mdp_count <- mdp_counts[mdp_row,mdp_col]
        background_factor <- background_table[corres_mitogene,mdp_col]
        corrected_mdp_count <- mdp_count - (background_factor*(corres_proportion)) * initial_mdp_count
        mdp_counts[mdp_row,mdp_col] <- corrected_mdp_count
      }
    }
  }
  for (mdp_col in 2:ncol(mdp_counts)) {
    mdp_counts[,mdp_col] <- as.integer(mdp_counts[,mdp_col])
  }
  return(mdp_counts)
}
