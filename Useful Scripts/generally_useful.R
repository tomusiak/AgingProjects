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

makeBAMS <- function(directory,paired_end) {
  fastq_list <- list.files(path=".",pattern="*.fastq",all.files=TRUE,full.names=FALSE)
  if (!paired_end) {
    for (fastq in 1:length(fastq_list)) {
      read_1 <- fastq_list[fastq]
      align(index="/home/atom/Desktop/Data/reference_index",readfile1=read_1,
            output_file=paste(c(str_sub(read_1,1,nchar(read_1)-6),".BAM"),collapse=""),
            nthreads = 4)
    }
  }
  if (paired_end) {
    list_of_indices <- rep(c(1,0),length(fastq_list)/2)
    for (fastq in 1:length(fastq_list)) {
      if(list_of_indices[fastq] == 1) {
        read_1 <- fastq_list[fastq]
        read_2 <- fastq_list[fastq+1]
        align(index="/home/atom/Desktop/Data/reference_index",readfile1=read_1,readfile2=read_2,
              output_file=paste(c(str_sub(read_1,1,nchar(read_1)-6),".BAM"),collapse=""),
              nthreads = 4)
      }
    }
  }
  bam_list <- list.files(path=".",pattern="*.BAM$",all.files=TRUE,full.names=FALSE)
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
  rownames(mitogene_count_matrix) <- sub('MT-', '', rownames(mitogene_count_matrix))
  write.csv(mitogene_count_matrix,"mitogene_count_matrix.csv") #Writes count matrix for easier future loading.
  return(mitogene_count_matrix)
}

getCountsMDP <- function(bam_list,paired_end) {
  gtffile <- file.path("/home/atom/Desktop/Data/mdp_atg_noATC.gtf")
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
  mdp_counts <- assays(se)$counts
  rownames(mdp_counts) <- sub('.Peptide', '', rownames(mdp_counts))
  mdp_counts <- mdp_counts[c(-593,-595),]
  keep <- rowSums((mdp_counts)) >= 200 #Removes genes with low counts.
  mdp_counts<- mdp_counts[keep,]
  write.csv(mdp_counts,"raw_counts_mdps.csv")
  return(mdp_counts)
}

deleteBAMs <- function() {
  bam_list <- list.files(path=".",pattern="*.BAM",all.files=TRUE,full.names=FALSE)
  for (bam in bam_list) {
    unlink(bam)
  }
}

deleteFASTQs <- function() {
  fastq_list <- list.files(path=".",pattern="*.fastq",all.files=TRUE,full.names=FALSE)
  for (fastq in fastq_list) {
    unlink(fastq)
  }
}

createCPGTable <- function(cpg_annotation, genome_annotation, upstream, downstream) {
  cpg_annotation <- cpg_annotation[,1:5]
  genome_annotation <- genome_annotation[genome_annotation$V3 == "gene" | genome_annotation$V3 == "ncRNA",]
  genome_annotation <- separate(genome_annotation, V9, c("gene_id","gene_type","gene_name"), ";")
  desired_columns <- c(1,3,4,5,7,11)
  genome_annotation <- genome_annotation[,desired_columns]
  genome_annotation <- separate(genome_annotation,gene_name,c("trash","trash2","gene_name")," ")
  genome_annotation <- genome_annotation[,c(1,2,3,4,5,8)]
  colnames(genome_annotation) <- c("chr","type","start","stop","strand","gene_name")
  colnames(cpg_annotation) <- c("chr","start","stop","strand","probe_name")
  cpg_table <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(cpg_table) <- c("cpg","gene")
  for (cpg in 1:nrow(cpg_annotation)) {
    cpg_chr <- cpg_annotation[cpg,1]
    cpg_start <- cpg_annotation[cpg,2]
    if(!is.na(cpg_start)) {
      relevant_genes <- genome_annotation[genome_annotation$chr == cpg_chr,]
      associated_genes <- relevant_genes[(relevant_genes$start - cpg_start <= downstream &
                                            relevant_genes$start - cpg_start >= 0) |
                                           (cpg_start - relevant_genes$start <= upstream &
                                              cpg_start - relevant_genes$start >= 0),]
      if (nrow(associated_genes) != 0) {
        associated_genes$cpg <- cpg_annotation[cpg,5]
        added_table <- associated_genes[,c(6,7)]
        cpg_table <- rbind(cpg_table,added_table)
      }
    }
  }
  return(cpg_table)
}

grabCPGs <- function(gene_name) {
  return(cpg_table[cpg_table$gene_name == gene_name,2])
}

getDiffMethylation <- function(gene_name,cpg_table) {
  gene <- beta_values[grabCPGs(gene_name),]
  gene <- data.frame(t(gene))
  gene$mean <- rowMeans(gene)
  gene$type <- all_data$type
  gene$type <- c(rep("Naive",7),rep("CM",7),rep("EM",7),rep("TEMRA",6))
  gene <- gene[,c(ncol(gene)-1,ncol(gene))]
  gene_summary <- getSummary(gene,"mean","type")
  return(gene_summary[,c(1,3)])
}

getDiffMethylationList <- function(gene_list, cpg_table) {
  diff_meth <- data.frame(matrix(ncol = 0, nrow = 4))
  diff_meth$type <- c("naive","central_memory","effector_memory","temra")
  for (gene in 1:length(gene_list)) {
    table <- getDiffMethylation(gene_list[gene],cpg_table)
    colnames(table) <- c("type",gene_list[gene])
    diff_meth <- cbind(diff_meth,table[,gene_list[gene],drop=FALSE])
  }
  return (diff_meth)
}

#Let's focus on these TFs.
getDiffMethylation2 <- function(gene_name,cpg_table) {
  grabCPGs <- function(gene_name) {
    return(cpg_table[cpg_table$gene_name == gene_name,2])
  }
  gene <- beta_values[grabCPGs(gene_name),]
  gene <- data.frame(t(gene))
  gene$type <- c(rep("Naive",7),rep("CM",7),rep("EM",7),rep("TEMRA",6))
  gene$name <- gene_name
  return (gene)
}

getDiffMethylationList2 <- function(gene_list, cpg_table) {
  diff_meth <- getDiffMethylation2(gene_list[1],cpg_table)
  for (gene in 2:length(gene_list)) {
    table <- getDiffMethylation2(gene_list[gene],cpg_table)
    diff_meth <- rbind.fill(diff_meth,table)
  }
  return (diff_meth)
}            nthreads = 4)
    }
  }
  if (paired_end) {
    list_of_indices <- rep(c(1,0),length(fastq_list)/2)
    for (fastq in 1:length(fastq_list)) {
      if(list_of_indices[fastq] == 1) {
        read_1 <- fastq_list[fastq]
        read_2 <- fastq_list[fastq+1]
        align(index="/home/atom/Desktop/Data/reference_index",readfile1=read_1,readfile2=read_2,
              output_file=paste(c(str_sub(read_1,1,nchar(read_1)-6),".BAM"),collapse=""),
              nthreads = 4)
      }
    }
  }
  bam_list <- list.files(path=".",pattern="*.BAM$",all.files=TRUE,full.names=FALSE)
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
  rownames(mitogene_count_matrix) <- sub('MT-', '', rownames(mitogene_count_matrix))
  write.csv(mitogene_count_matrix,"mitogene_count_matrix.csv") #Writes count matrix for easier future loading.
  return(mitogene_count_matrix)
}

getCountsMDP <- function(bam_list,paired_end) {
  gtffile <- file.path("/home/atom/Desktop/Data/mdp_atg_noATC.gtf")
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
  mdp_counts <- assays(se)$counts
  rownames(mdp_counts) <- sub('.Peptide', '', rownames(mdp_counts))
  mdp_counts <- mdp_counts[c(-593,-595),]
  keep <- rowSums((mdp_counts)) >= 200 #Removes genes with low counts.
  mdp_counts<- mdp_counts[keep,]
  write.csv(mdp_counts,"raw_counts_mdps.csv")
  return(mdp_counts)
}

deleteBAMs <- function() {
  bam_list <- list.files(path=".",pattern="*.BAM",all.files=TRUE,full.names=FALSE)
  for (bam in bam_list) {
    unlink(bam)
  }
}

deleteFASTQs <- function() {
  fastq_list <- list.files(path=".",pattern="*.fastq",all.files=TRUE,full.names=FALSE)
  for (fastq in fastq_list) {
    unlink(fastq)
  }
}

createCPGTable <- function(cpg_annotation, genome_annotation, upstream, downstream) {
  cpg_annotation <- cpg_annotation[,1:5]
  genome_annotation <- genome_annotation[genome_annotation$V3 == "gene" | genome_annotation$V3 == "ncRNA",]
  genome_annotation <- separate(genome_annotation, V9, c("gene_id","gene_type","gene_name"), ";")
  desired_columns <- c(1,3,4,5,7,11)
  genome_annotation <- genome_annotation[,desired_columns]
  genome_annotation <- separate(genome_annotation,gene_name,c("trash","trash2","gene_name")," ")
  genome_annotation <- genome_annotation[,c(1,2,3,4,5,8)]
  colnames(genome_annotation) <- c("chr","type","start","stop","strand","gene_name")
  colnames(cpg_annotation) <- c("chr","start","stop","strand","probe_name")
  cpg_table <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(cpg_table) <- c("cpg","gene")
  for (cpg in 1:nrow(cpg_annotation)) {
    cpg_chr <- cpg_annotation[cpg,1]
    cpg_start <- cpg_annotation[cpg,2]
    if(!is.na(cpg_start)) {
      relevant_genes <- genome_annotation[genome_annotation$chr == cpg_chr,]
      associated_genes <- relevant_genes[(relevant_genes$start - cpg_start <= downstream &
                                           relevant_genes$start - cpg_start >= 0) |
                                           (cpg_start - relevant_genes$start <= upstream &
                                              cpg_start - relevant_genes$start >= 0),]
      if (nrow(associated_genes) != 0) {
        associated_genes$cpg <- cpg_annotation[cpg,5]
        added_table <- associated_genes[,c(6,7)]
        cpg_table <- rbind(cpg_table,added_table)
      }
      }
    }
  return(cpg_table)
}

getDiffMethylation <- function(gene_name,cpg_table) {
  grabCPGs <- function(gene_name) {
    return(cpg_table[cpg_table$gene_name == gene_name,2])
  }
  gene <- beta_values[grabCPGs(gene_name),]
  gene <- data.frame(t(gene))
  gene$mean <- rowMeans(gene)
  gene$type <- all_data$type
  gene$type <- c(rep("Naive",7),rep("CM",7),rep("EM",7),rep("TEMRA",6))
  gene <- gene[,c(ncol(gene)-1,ncol(gene))]
  gene_summary <- getSummary(gene,"mean","type")
  return(gene_summary[,c(1,3)])
}

getDiffMethylationList <- function(gene_list, cpg_table) {
  diff_meth <- data.frame(matrix(ncol = 0, nrow = 4))
  diff_meth$type <- c("naive","central_memory","effector_memory","temra")
  for (gene in 1:length(gene_list)) {
    table <- getDiffMethylation(gene_list[gene],cpg_table)
    colnames(table) <- c("type",gene_list[gene])
    diff_meth <- cbind(diff_meth,table[,gene_list[gene],drop=FALSE])
  }
  return (diff_meth)
}

#Let's focus on these TFs.
getDiffMethylation2 <- function(gene_name,cpg_table) {
  grabCPGs <- function(gene_name) {
    return(cpg_table[cpg_table$gene_name == gene_name,2])
  }
  gene <- beta_values[grabCPGs(gene_name),]
  gene <- data.frame(t(gene))
  gene$type <- c(rep("Naive",7),rep("CM",7),rep("EM",7),rep("TEMRA",6))
  gene$name <- gene_name
  return (gene)
}

getDiffMethylationList2 <- function(gene_list, cpg_table) {
  diff_meth <- getDiffMethylation2(gene_list[1],cpg_table)
  for (gene in 2:length(gene_list)) {
    table <- getDiffMethylation2(gene_list[gene],cpg_table)
    diff_meth <- rbind.fill(diff_meth,table)
  }
  return (diff_meth)
}
