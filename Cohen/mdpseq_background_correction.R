#Alan Tomusiak

#Contains functions useful for MDPSeq background correction.
# Functions by processing changes relative to the mean in the mitochondrial genes (hereafter,
# mitogenes), parsing an "encompass table" that maps MDPs to mitogenes and their proportion of
# overlap, and then subtracting the signal in the mitogenes.
# Can be useful in datasets where large changes in mitogene expression can skew MDPSeqs, and
# intuitively the code has seemed to perform well in practice.

#A few caveats:
# 1) Does not incorporate strand-specific information. In fact, it is likely that this will
# provide false signal with strand-specific data. If desired, the code could actually be fairly easily
# altered to use the strand-specific data effectively.
#
# 2) If an MDP overlaps two mitogenes, it will only background correct for the larger overlap.
#
# 3) Due to the small size of tRNAs, background correction will struggle to correctly 
# account for changes in tRNA expression. 

#Loading in necessary packages.
library(tidyr)
library(tidyverse)

#Adjustable parameter that will assist in estimation of how much an MDP signal will
# impact mitogene reads. Particularly important if attempting to background correct
# tRNAs.
read_length <- 150

#Function to create an "encompass table" mapping MDPs to mitogenes and calculating
# percentage of overlap.
generateEncompassTable <- function(gtf_mdp, gtf_mitochondrial) {
  #Pre-processing MDP gtf. Removes unnecessary information.
  mdp_gtf <- gtf_mdp
  colnames(mdp_gtf) <-
    c("chr","idc","idc2", "start","end",
      "idc3","sense","idc4","gene")
  mdp_gtf <- separate(mdp_gtf, gene, c("gene", "transcript"), ";")
  mdp_gtf <- separate(mdp_gtf, gene, c("trash", "mdp_id"), ">Peptide")
  mdp_db <- mdp_gtf[, c("start", "end", "sense", "mdp_id")]
  mdp_db <- mdp_db[-c(62, 425), ]
  
  #Pre-processing mitochondrial GTF. Removes unnecessary information.
  mitogene_gtf <- gtf_mitochondrial
  mitogene_db <- mitogene_gtf[, c("V3", "V4", "V5", "V7", "V9")]
  mitogene_db <- mitogene_db[mitogene_db$V3 == "exon", ]
  colnames(mitogene_db) <-
    c("exon", "start", "end", "sense", "mitogene_name")
  mitogene_db$mitogene_name <-
    c(
      "TF","RNR1","TV","RNR2","TL1","ND1","TI",
      "TQ","TM","ND2","TW","TA","TN","TC","TY",
      "CO1","TS","TD","CO2","TK","ATP8","ATP6",
      "CO3","TG","ND3","TR","ND4L","ND4","TH",
      "TS2","TL2","ND5","ND6","TE","CYB","TT","TP"
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
  # The way this function works is that it effectively checks for the four overlap possibilities -
  # the MDP is entirely in the mitogene, the MDP is entirely out of the mitogene, the MDP is
  # partially inside the mitogene (5' end), or it is partially inside the mitogene (3' end).
  # The "proportion" parameter will be useful later for calculating the amount of signal that
  # the MDP is providing to the mitogene reads.
  for (mdp_row in 1:nrow(mdp_db)) { # Loops through all MDPs.
    mdp_s <- mdp_db[mdp_row, "start"]
    mdp_e <- mdp_db[mdp_row, "end"]
    mdp_size <- mdp_e-mdp_s # MDP size useful for calculating % overlap.
    for (mito_row in 1:nrow(mitogene_db)) { # Loops through all mitogenes.
      mito_s <- mitogene_db[mito_row, "start"]
      mito_e <- mitogene_db[mito_row, "end"]
      mito_size <- mito_e-mito_s # Mitogene size useful for calculating % overlap.
      if ((mdp_s < mito_s & # No overlap.
           mdp_s < mito_e & mdp_e < mito_s & mdp_e < mito_e) |
          (mdp_s > mito_s &
           mdp_s > mito_e & mdp_e > mito_s & mdp_e > mito_e)) {
        next;
      }
      #Checks for partial overlap on 5' end.
      if (mdp_s < mito_s &
          mdp_s < mito_e & mdp_e > mito_s & mdp_e < mito_e) {
        overlap <- (mdp_e - mito_s) / (mdp_e - mdp_s)
        encompass_table <- encompass_table %>% add_row(mdp = mdp_db[mdp_row, "mdp_id"],
                                    mitogene = mitogene_db[mito_row, "mitogene_name"],
                                    perc_overlap = overlap, 
                                    proportion=overlap*(((mdp_size)+read_length)/mito_size))
      } 
      #Checks for partial overlap on 3' end.
      if (mdp_s > mito_s &
          mdp_s < mito_e & mdp_e > mito_s & mdp_e > mito_e) {
        overlap <- (mito_e - mdp_s) / (mdp_e - mdp_s)
        encompass_table <- encompass_table %>% add_row(mdp = mdp_db[mdp_row, "mdp_id"],
                                    mitogene = mitogene_db[mito_row, "mitogene_name"],
                                    perc_overlap = overlap,
                                    proportion=overlap*(((mdp_size)+read_length)/mito_size))
      } 
      #Checks for 100% overlap.
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
  
  #A somewhat hacky solution to deal with the case where an MDP overlaps multiple mitogenes.
  # This code effectively creates a new "encompass table" that only takes in the mitogene with
  # the highest amount of overlap.
  encompass_table<-encompass_table[!is.na(encompass_table$mdp),]
  encompass_table_new <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(encompass_table_new) <- c("mdp", "mitogene", "perc_overlap","proportion")
  encompass_table_new$mdp <- as.character(encompass_table_new$mdp)
  encompass_table_new$mitogene <- as.character(encompass_table_new$mitogene)
  encompass_table_new$perc_overlap <-
    as.double(encompass_table_new$perc_overlap)
  encompass_table_new$proportion <-
    as.double(encompass_table_new$proportion)
  for (mdp_row in 1:nrow(mdp_db)) {
    mdp <- mdp_db[mdp_row,"mdp_id"]
    encompassing_mdps <- encompass_table[encompass_table$mdp == mdp,]
    if (nrow(encompassing_mdps)!=0) {
      encompassing_mdps <- encompassing_mdps[order(encompassing_mdps$proportion,decreasing=FALSE),]
      first_mdp <- encompassing_mdps[1,]
      encompass_table_new <- rbind(encompass_table_new,first_mdp)
    }
  }
  encompass_table_new <- drop_na(encompass_table_new)
  
  #Spitting out the data.
  return (encompass_table_new)
}

#A fairly straightforward function that determines the amount of background signal in each
# mitogene by determining deviation from the median. This will lead to a small transformation of
# the dataset downstream, particularly if there are massive outliers.
determineBackgroundSignal <- function(mito_counts) {
  median <- apply(mito_counts, 1, median, na.rm=T)
  mito_counts$median <- median
  for (mito_column in 1:(ncol(mito_counts)-1)) {
      mito_counts[,mito_column] <- (mito_counts[,mito_column] - mito_counts[,ncol(mito_counts)]) /
        ((mito_counts[,ncol(mito_counts)]))
  }
  mito_counts <- mito_counts[,-ncol(mito_counts)]
  return (mito_counts)
}

#Performs background correction of the dataset by subtracting signal generated by
#mitogenes from the signals of their corresponding MDPs.
performBackgroundCorrection <- function(background_table,mdp_counts,encompass_table) {
  mdp_names <- rownames(mdp_counts)
  for (mdp_row in 1:nrow(mdp_counts)) { #Loops through MDPs.
    #Calculates proportion information for the corresponding mitogene of each MDP.
    median_mdp <- apply(mdp_counts[mdp_row,], 1, median, na.rm=T)
    mdp <- mdp_names[mdp_row]
    corres_mitogene <- encompass_table$mitogene[encompass_table$mdp==mdp]
    corres_overlap <- encompass_table$perc_overlap[encompass_table$mdp==mdp]
    corres_proportion <- 1-(encompass_table$proportion[encompass_table$mdp==mdp])
    if(length(corres_proportion) != 0) {
      #If MDP signal on the mitogene is overwhelming, no background
      # correction is performed.
      if (corres_proportion < 0) { 
        corres_proportion <- 0
      } else {
        for (mdp_col in 1:length(mdp_counts[mdp_row,])) {
          mdp_count <- mdp_counts[mdp_row,mdp_col]
          background_factor <- background_table[corres_mitogene,mdp_col]
          if (!is.nan(background_factor)) {
            # Subtracts out mitogene signal from MDP signal.
            corrected_mdp_count <- mdp_count - (background_factor*(corres_proportion)) * median_mdp
            mdp_counts[mdp_row,mdp_col] <- corrected_mdp_count
            if(mdp_counts[mdp_row,mdp_col] < 0) {
              mdp_counts[mdp_row,mdp_col] <- 0
            }
          }
        }
      }
    }
  }
  for (mdp_col in 1:ncol(mdp_counts)) {
    mdp_counts[,mdp_col] <- as.integer(mdp_counts[,mdp_col])
  }
  return(mdp_counts)
}
