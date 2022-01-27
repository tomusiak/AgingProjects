library(stringr)

#ripped this from someone on the internet
getSummary <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                       conf.interval=.95, .drop=TRUE, removeOutliers=F, logt=F) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         
                         if (removeOutliers) {
                           print('pre')
                           print(xx[[col]])
                           #xx[[col]] <- boxB(xx[[col]])
                           print('post')
                           #print(exp(rm.outlier(log(xx[[col]]), fill = FALSE, median = FALSE, opposite = FALSE)))
                           print(scores(log(xx[[col]]), type="t", prob=0.90) )
                           #print(xx[[col]][boxB(xx[[col]], method="resistant", logt=logt, k=2.5)[["outliers"]]])
                         }
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

makeBAMS <- function(directory) {
  fastq_list <- list.files(path=".",pattern="*.fastq",all.files=TRUE,full.names=FALSE)
  for (fastq in fastq_list) {
    align(index="/home/atom/Desktop/Data/reference_index",readfile1=file,
          output_file=paste(c(str_sub(fastq,1,nchar(fastq)-6),".BAM"),collapse=""),
          nthreads = 4)
    unlink(fastq)
  }
  bam_list <- list.files(path=".",pattern="*.BAM",all.files=TRUE,full.names=FALSE)
  return (bam_list)
}

getCountsMitochondrial <- function(bam_list,paired_end) {
  mitogene_counts <- featureCounts(bam_list,annot.ext="~/Desktop/Data/mitochondria_db.gtf",
                                   isGTFAnnotationFile=TRUE,nthreads=5) #Count matrix generation
  mitogene_count_matrix <- mitogene_counts$counts
  rownames(mitogene_count_matrix) <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1",
                                       "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2",
                                       "MT-TW", "MT-TA", "MT-TN", "MT-TC","MT-TY",
                                       "MT-CO1", "MT-TS","MT-TD", "MT-CO2","MT-TK",
                                       "MT-ATP8","MT-ATP6","MT-CO3", "MT-TG", "MT-ND3",
                                       "MT-TR","MT-ND4L","MT-ND4","MT-TH","MT-TS2",
                                       "MT-TL2","MT-ND5","MT-ND6","MT-TE","MT-CYB",
                                       "MT-TT","MT-TP")
  write.csv(mitogene_count_matrix,"mitogene_count_matrix.csv") #Writes count matrix for easier future loading.
  for (file in bam_list) {
    unlink(file)
  }
  return(mitogene_count_matrix)
}

getCountsMDP <- function(bam_list,paired_end) {
  gtffile <- file.path(getwd(),"/home/atom/Desktop/Data/mdp_atg_noATC.gtf")
  txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
  mitogene_gtf <-
    read.delim(
      "/home/atom/Desktop/Data/mitochondria_db.gtf",
      header = FALSE,
      comment.char = "#"
    )
  ebg <- exonsBy(txdb, by="gene")
  se <- summarizeOverlaps(features=ebg, reads=bam_list,
                          mode="Union",
                          singleEnd=!paired_end,
                          ignore.strand=FALSE,
                          fragments=FALSE,
                          inter.feature=FALSE) ###this will count each)
  saveRDS(se, "mdp_se.rds")
  mdp_counts <- assays(se_cancer)$counts
  rownames(mdp_counts) <- sub('.Peptide', '', rownames(mdp_counts))
  mdp_counts <- mdp_counts[c(-593,-595),]
  keep <- rowSums((mdp_counts)) >= 200 #Removes genes with low counts.
  mdp_counts<- mdp_counts[keep,]
  write.csv(mdp_counts,"raw_counts_mdps.csv")
  for (file in bam_list) {
    unlink(file)
  }
  return(mdp_counts)
}