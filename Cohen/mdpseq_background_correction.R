#load packages
library(tidyr)
library(tidyverse)

read_length <- 150

generateEncompassTable <- function(gtf_mdp, gtf_mitochondrial) {
  mdp_gtf <- gtf_mdp
  mitogene_gtf <- gtf_mitochondrial
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
  mdp_gtf <- separate(mdp_gtf, gene, c("trash", "mdp_id"), ">Peptide")
  mdp_db <- mdp_gtf[, c("start", "end", "sense", "mdp_id")]
  mdp_db <- mdp_db[-c(62, 425), ]
  
  #Pre-processing mitochondrial GTF.
  mitogene_db <- mitogene_gtf[, c("V3", "V4", "V5", "V7", "V9")]
  mitogene_db <- mitogene_db[mitogene_db$V3 == "exon", ]
  colnames(mitogene_db) <-
    c("exon", "start", "end", "sense", "mitogene_name")
  mitogene_db$mitogene_name <-
    c(
      "TF",
      "RNR1",
      "TV",
      "RNR2",
      "TL1",
      "ND1",
      "TI",
      "TQ",
      "TM",
      "ND2",
      "TW",
      "TA",
      "TN",
      "TC",
      "TY",
      "CO1",
      "TS",
      "TD",
      "CO2",
      "TK",
      "ATP8",
      "ATP6",
      "CO3",
      "TG",
      "ND3",
      "TR",
      "ND4L",
      "ND4",
      "TH",
      "TS2",
      "TL2",
      "ND5",
      "ND6",
      "TE",
      "CYB",
      "TT",
      "TP"
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
                                    perc_overlap = overlap, proportion=overlap*(((mdp_size)+read_length)/mito_size))
      } 
      if (mdp_s > mito_s &
          mdp_s < mito_e & mdp_e > mito_s & mdp_e > mito_e) {
        overlap <- (mito_e - mdp_s) / (mdp_e - mdp_s)
        encompass_table <- encompass_table %>% add_row(mdp = mdp_db[mdp_row, "mdp_id"],
                                    mitogene = mitogene_db[mito_row, "mitogene_name"],
                                    perc_overlap = overlap,proportion=overlap*(((mdp_size)+read_length)/mito_size))
      } 
      if (mdp_s > mito_s &
          mdp_s < mito_e & mdp_e > mito_s & mdp_e < mito_e) {
        overlap <- 1
        encompass_table <- encompass_table %>%
          add_row(mdp = mdp_db[mdp_row, "mdp_id"],
                  mitogene = mitogene_db[mito_row, "mitogene_name"],
                  perc_overlap = overlap,
                  proportion = overlap*((mdp_size)+read_length)/mito_size)
      } 
    }
  }
  encompass_table<-encompass_table[!is.na(encompass_table$mdp),]
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
            mdp_row_2 <- 1
            mdp_row <- 1
          }
          if (outer_overlap == inner_overlap) {
            encompass_table <- encompass_table[-mdp_row_2,]
            mdp_row_2 <- 1
            mdp_row <- 1
          }
          if (outer_overlap < inner_overlap) {
            encompass_table <- encompass_table[-mdp_row,]
            mdp_row_2 <- 1
            mdp_row <- 1
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
        ((mito_counts[,ncol(mito_counts)]))
  }
  mito_counts <- mito_counts[,-ncol(mito_counts)]
  return (mito_counts)
}

performBackgroundCorrection <- function(background_table,mdp_counts,encompass_table) {
  mdp_names <- rownames(mdp_counts)
  for (mdp_row in 1:nrow(mdp_counts)) {
    average_mdp <- mean(mdp_counts[mdp_row,])
    mdp <- mdp_names[mdp_row]
    corres_mitogene <- encompass_table$mitogene[encompass_table$mdp==mdp]
    corres_overlap <- encompass_table$perc_overlap[encompass_table$mdp==mdp]
    corres_proportion <- 1-(encompass_table$proportion[encompass_table$mdp==mdp])
    if (length(corres_proportion) != 0) {
      if (corres_proportion < 0) {
        corres_proportion <- 0
      } 
    }
    if(length(corres_mitogene) != 0) {
      for (mdp_col in 1:length(mdp_counts[mdp_row,])) {
        mdp_count <- mdp_counts[mdp_row,mdp_col]
        background_factor <- background_table[corres_mitogene,mdp_col]
        corrected_mdp_count <- mdp_count - (background_factor*(corres_proportion)) * average_mdp
        mdp_counts[mdp_row,mdp_col] <- corrected_mdp_count
        if(mdp_counts[mdp_row,mdp_col] < 0) {
          mdp_counts[mdp_row,mdp_col] <- 0
        }
      }
    }
  }
  for (mdp_col in 1:ncol(mdp_counts)) {
    mdp_counts[,mdp_col] <- as.integer(mdp_counts[,mdp_col])
  }
  return(mdp_counts)
}
