#--------------------------------------------------------#
#-- Frances' longitudinal data from the Lostine/Asotin --#
#--------------------------------------------------------#

#-- Read in and prep dataframe --#

data <- read.csv(file.choose(), header = T, sep = "\t")
#-- read in the csv version of the longitudinal sampling data. The file.choose part of the command below opens a navigation window that --#
#-- will let you navigate directly to the dataset --#
head(data) 
#-- Shows the first five rows of the dataset. I use this to be sure that the dataset was read in with the correct delimiting.
data$WADDL.Movi.PCR <- ifelse(data$WADDL.Movi.PCR == "Indeterminant", "Indeterminate", as.character(data$WADDL.Movi.PCR))
data$WADDL.Movi.PCR <- factor(data$WADDL.Movi.PCR)
#-- These two lines correct the indeterminate/indeterminant inconsistency in the dataset. 
#-- After these lines, all indeterminate shedders are labeled "Indeterminate".
table(data$WADDL.Movi.PCR)

#-- build and explore determination cohort --#

Animal.tab <- table(data$Animal_ID)
#-- Tables number of observations for each Animal ID and stores table in "Animal.tab" object
AnimalByOutcome.tab <- table(data$Animal_ID, data$WADDL.Movi.PCR)
Animal.tab[which(Animal.tab >= 4)]
#-- Extracts Animal IDs with 4 or more samples
Animal.tab[which(Animal.tab == 3)]
#-- Extracts Animal IDs with 3 samples

ManySamplesOutcomes.tab <- AnimalByOutcome.tab[which(Animal.tab >=4), ]
ManySamplesProps <- ManySamplesOutcomes.tab[ ,1] / (ManySamplesOutcomes.tab[ ,1] + ManySamplesOutcomes.tab[ ,4])
ManyAnimalSamples.full <- subset(data, Animal_ID %in% names(Animal.tab[which(Animal.tab >= 4)]))

#--------------------------------------------------------------------------------#
#-- Histogram of proportion of positive tests for sheep with 4 or more samples --#
#--------------------------------------------------------------------------------#

hist(ManySamplesProps, breaks = 10, col = "grey80", xlab = "Proportion of PCR-positive tests for \n the 12 animals tested 4 or more times", main = "")
#-- build a histogram of the proportion of tests that were positive for animals tested 4 or more times

#-- classify animals in determination cohort into carriage states ----#
#-- place cut line at .6: anybody above .6 is classified as chronic --#
determination.cohort <- subset(data, Animal_ID %in% names(which(Animal.tab >= 4)))
Chronic.orig <- c("01LO34", "04LO12", "99L09")
Intermittent.orig <- c("04LO73", "05LO68", "10LO33", "12LO56", "12LO75")
determination.cohort$CarriageState <- ifelse(determination.cohort$Animal_ID %in% Chronic.orig, "chronic", ifelse(determination.cohort$Animal_ID %in% Intermittent.orig, "inter", "not.chronic"))
ManyAnimalSamples.full$CarriageState <- ifelse(ManyAnimalSamples.full$Animal_ID %in% Chronic.orig, "chronic", ifelse(ManyAnimalSamples.full$Animal_ID %in% Intermittent.orig, "inter", "not.chronic"))


#---------------------------------#
#-- JAGS mixture trials ----------#
#---------------------------------#
inits = function() {
    lambda1 = mean(eyesdata$y[1:30]) +rnorm(1,0,.01)
    theta = mean(eyesdata$y[31:48]) - lambda1
    sigma2 = var(eyesdata$y[1:30])
    return(list(lambda = c(lambda1,NA), theta = theta, tau = 1/sigma2, pi  = c(30, 48-30)/48))
  }

#--------------------------------------------#
#-- R mixture discriminant analysis trials --#
#--------------------------------------------#
require(mda)
mda.test <- mda(CarriageState ~ WADDL.Movi.PCR + Movipneumonia.ELISA, data = ManyAnimalSamples.full)
mda.preds <- predict(mda.test)

new.data <- subset(data, select = c(WADDL.Movi.PCR, Movipneumonia.ELISA))
new.preds <- predict(mda.test, new.data, type = "posterior")

require(MASS)
lda.test <- lda(CarriageState ~ WADDL.Movi.PCR + Movipneumonia.ELISA, data = ManyAnimalSamples.full)
lda.test$svd
lda.preds <- predict(lda.test)$class

qda.test <- qda(CarriageState ~ WADDL.Movi.PCR + Movipneumonia.ELISA, data = ManyAnimalSamples.full)
qda.test$svd
qda.preds <- predict(qda.test)$class

x <- cbind(ManyAnimalSamples.full$CarriageState, as.character(lda.preds), as.character(qda.preds), as.character(mda.preds))
model.sens <- 15/16
model.spec <- 32/39

ThreeSamplesOutcomes.tab <- AnimalByOutcome.tab[which(Animal.tab == 3), ]
Chronic.3 <- c("01LO32", "04LO74", "13LO79")
Intermittent.3 <- c("03LO50", "04LO65", "11AS07", "12AS03")
Neg.3 <- c("04LO49", "05LO53", "06LO61", "12AS16")

test.cohort <- subset(data, Animal_ID %in% names(which(Animal.tab == 3)))
test.cohort$CarriageState <- ifelse(test.cohort$Animal_ID %in% Chronic.3, "chronic", ifelse(test.cohort$Animal_ID %in% Intermittent.3, "inter", "neg"))
