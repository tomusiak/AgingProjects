#This program creates a new epigenetic clock by training on CpGs that were identified to be not associated with
# T cell differentiation in a previous analysis. The clock is based on elastic net.

#Some useful code.
setwd("Data/") #Sets directory.

#Reading in necessary libraries.
library(glmnet)
library(readr)
library(tidyr)
library(caret)
library(MASS)
library(glmnetUtils)
library(random)
library(Rfast)
library(dplyr)
library(limma)
library(MASS)
library(sva)

#Reading in database data for training and testing.
ml_cpg_table <- data.table::fread("ClockConstruction/healthy_cpg_table.csv",header=TRUE) %>% 
  as.data.frame()
row.names(ml_cpg_table) <- ml_cpg_table$V1
ml_cpg_table <- ml_cpg_table[,-1]
permitted_cpgs <- read.csv("ClockConstruction/nodiff_cpgs.csv",row.names=1)[,1]
ml_sample_table <- read.csv("ClockConstruction/healthy_sample_table.csv",row.names=1)
# ml_cpg_table <- ml_cpg_table[rownames(ml_cpg_table) %in% permitted_cpgs,]

ml_cpg_table[ml_cpg_table == "NULL"] <- NA
ml_cpg_table[ml_cpg_table == "null"] <- NA
ml_cpg_table <- ml_cpg_table[rowSums(is.na(ml_cpg_table))==0,]
ml_cpg_table <- drop_na(ml_cpg_table)
ml_sample_table <- ml_sample_table[ml_sample_table$ID %in% colnames(ml_cpg_table),]

cpgs <- rownames(ml_cpg_table)
samples <- colnames(ml_cpg_table)
ml_cpg_table <- as.matrix(sapply(ml_cpg_table, as.numeric)) 
ml_cpg_table <- ComBat(ml_cpg_table,ml_sample_table$Author,mean.only=TRUE)

#ml_cpg_table <- scale(ml_cpg_table)

#Annotating and splitting data.
ml_cpg_table_rotated <- data.frame(Rfast::transpose(ml_cpg_table))
colnames(ml_cpg_table_rotated) <- cpgs
rownames(ml_cpg_table_rotated) <- samples

ml_cpg_table_rotated <- ml_cpg_table_rotated[ rownames(ml_cpg_table_rotated)%in% ml_sample_table$ID ,]
ml_cpg_table_rotated$Age <- ml_sample_table$Age

# adult_age <- 20
# 
# ml_cpg_table_rotated[ml_cpg_table_rotated$Age > adult_age]$Age <- ((ml_cpg_table_rotated$Age)-adult_age)/(adult_age+1)
# ml_cpg_table_rotated[ml_cpg_table_rotated$Age <= adult_age]$Age <- log(ml_cpg_table_rotated$Age+1)-log(adult_age+1)


#Reading in validation table of data from Jonkman 2022. Will be useful for comparing clock performance.
validation_cpg_table <- data.table::fread("ClockConstruction/validation_cpg_table.csv",
                                           header=TRUE) %>%
                                           as.data.frame()
row.names(validation_cpg_table) <- validation_cpg_table$V1
validation_cpg_table <- validation_cpg_table[,-1]
validation_cpg_table <- validation_cpg_table[rownames(validation_cpg_table) %in%
                                               colnames(ml_cpg_table_rotated),]
validation_cpg_table_rotated <- data.frame(t(validation_cpg_table))
validation_cpg_table_rotated <- validation_cpg_table_rotated[,colnames(validation_cpg_table_rotated) %in%
                                                               permitted_cpgs]

#Manually separating out datasets into training and test datasets.
all_authors <- (ml_sample_table %>% group_by(Author) %>% dplyr::summarize(Authors=n()))$Author
training_set_chosen <- sample(c(1:length(all_authors)), (length(all_authors)*2/3))
training_authors <- all_authors[training_set_chosen]
test_authors <- all_authors[-training_set_chosen]
training_samples <- ml_sample_table[ml_sample_table$Author %in% training_authors,]
training_samples <- training_samples$ID
test_samples <- ml_sample_table[ml_sample_table$Author %in% test_authors,]
test_samples <- test_samples$ID
training_set <- ml_cpg_table_rotated[training_samples,]
test_set <- ml_cpg_table_rotated[test_samples,]

#delete this after..
training_set_alternate_slate <- sample(c(1:length(rownames(ml_cpg_table_rotated))), (length(rownames(ml_cpg_table_rotated))*2/3))
training_set_alternate <- ml_cpg_table_rotated[training_set_alternate_slate,]
test_set_alternate <- ml_cpg_table_rotated[-training_set_alternate_slate,]
training_set_alternate <- na.omit(training_set)
test_set_alternate <- na.omit(test_set_alternate)

#Training the model.

fit_control <- trainControl(
      method = "repeatedcv",
      number = 2,
      repeats = 2,
      verboseIter=TRUE)

model <-train(Age ~ .,
      data = training_set_alternate,
      method = "glmnet",
      tuneLength = 25,
      preProcess= c("center","scale"),
      trControl = fit_control)

predicted_age <- predict(model,test_set_alternate)
predicted_age <- data.frame(predicted_age)
predicted_age$ID <- rownames(predicted_age)
sample_table_2 <- merge(ml_sample_table,predicted_age,by="ID")

#Plotting accuracy of the model.
fit = lm(predicted_age ~ Age, sample_table_2)
ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
  geom_point(aes(color=sample_table_2$Author)) +
  stat_smooth(method = "lm", col = "red") +
  theme_classic() +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "ME = ",mean(sqrt(fit$residuals^2)),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
# 
# 
# #Identifying which CpGs are important
# important_cpgs<-varImp(model)$importance
# important_cpgs <- important_cpgs[order(important_cpgs$Overall,decreasing=TRUE),,drop=FALSE]
# 
# 
# #Fitting the Jonkman useful datasets to the model and exporting them to analyze whether they work for
# # differentiation.
# fitted_validation <- predict(model,validation_cpg_table_rotated)
# write.csv(fitted_validation,"ClockConstruction/predictions.csv")

training_set <- na.omit(training_set)
test_set <- na.omit(test_set)
validation_set <- na.omit(validation_cpg_table_rotated)

write.csv(training_set_alternate,"ClockConstruction/training_set.csv")
write.csv(test_set_alternate,"ClockConstruction/test_set.csv")
write.csv(validation_set, "ClockConstruction/validation_set.csv")


