source("AgingProjects/Useful Scripts/generally_useful.R")
setwd("Data/") #Sets directory.

library(readr)
library(tidyr)
library(caret)
library(glmnet)
library(MASS)


sample_table <- read.csv("ClockConstruction/sample_table.csv",row.names=1)
cpg_table <- read.csv("ClockConstruction/cpg_table.csv",row.names=1)
cpg_table_rotated <- data.frame(t(cpg_table))


validation_sample_table <- read.csv("ClockConstruction/validation_sample_table.csv",row.names=1)
validation_cpg_table <- read.csv("ClockConstruction/validation_cpg_table.csv",row.names=1)
validation_cpg_table <- validation_cpg_table[rownames(validation_cpg_table) %in% 
                                               rownames(cpg_table),]
validation_cpg_table_rotated <- data.frame(t(validation_cpg_table))
validation_cpg_table_rotated <- validation_cpg_table_rotated[,colnames(cpg_table_rotated[,1:227605])]
small_validation <- validation_cpg_table_rotated[,220000:227605]

train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              search = "random",
                              verboseIter = TRUE)
cpg_table_rotated$Age <- sample_table$Age
small_test <- cpg_table_rotated[,220000:227606]
training <-train(Age ~ .,
      data = small_test,
      method = "glmnet",
      preProcess = c("center", "scale"),
      tuneLength = 25,
      trControl = train_control)
fitted <- predict(training)
sample_table_2 <- cbind(sample_table,fitted)
p <- ggplot(sample_table_2,aes(x=Age,y=fitted)) +
  geom_point() + theme_classic() + geom_smooth(method="lm")

lm_eqn <- function(df){
  m <- lm(Age ~ fitted, df);
  eq <- substitute(italic(Age) == a + b %.% italic(fitted)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(sample_table_2), parse = TRUE)

ggplotRegression <- function(fit){
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    theme_classic() +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(lm(Age ~ fitted, sample_table_2))

#curious
fitted_validation <- predict(training,small_validation)
write.csv(fitted_validation,"ClockConstruction/predictions.csv")
