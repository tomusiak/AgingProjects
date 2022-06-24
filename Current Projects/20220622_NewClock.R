
source("AgingProjects/Useful Scripts/generally_useful.R")
setwd("Data/") #Sets directory.

library(readr)
library(tidyr)
garaud_cpgs <- data.frame(read_table2("Garaud2017/GSE71825_series_matrix.txt", 
                                      col_names = FALSE, comment = "!"))
colnames(garaud_cpgs) <- garaud_cpgs[1,]
rownames(garaud_cpgs) <- garaud_cpgs[,1]
garaud_cpgs <- garaud_cpgs[-1,-1]
garaud_cpgs <- as.numeric(garaud_cpgs)
garaud_cpgs <- drop_na(garaud_cpgs)
garaud_cpgs <- as.character(garaud_cpgs)
garaud_cpgs<- as.numeric(garaud_cpgs)
