summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
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