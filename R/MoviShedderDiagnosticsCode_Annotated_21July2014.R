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

#--------------------------------------------------------------------------------#
#-- Histogram of proportion of positive tests for sheep with 4 or more samples --#
#--------------------------------------------------------------------------------#

hist(ManySamplesProps, breaks = 10, col = "grey80", xlab = "Proportion of PCR-positive tests for \n the 12 animals tested 4 or more times", main = "")
  #-- build a histogram of the proportion of tests that were positive for animals tested 4 or more times

#-- classify animals in determination cohort into carriage states ----#
#-- place cut line at .6: anybody above .6 is classified as chronic --#
determination.cohort <- subset(data, Animal_ID %in% names(which(Animal.tab >= 4)))
Chronic <- c("01LO34", "04LO12", "99L09")
determination.cohort$CarriageState <- ifelse(determination.cohort$Animal_ID %in% Chronic, "chronic", "not.chronic")

PropPos <- matrix(NA, ncol = 4, nrow = length(levels(factor(determination.cohort$Animal_ID))))
TrueStatus <- rep(NA, length(levels(factor(determination.cohort$Animal_ID))))
for(i in 1:4){
  tests.to.use <- subset(determination.cohort, TestNumber <= i)
  for(j in 1:length(levels(factor(determination.cohort$Animal_ID)))){
    anim.spec.data <- subset(tests.to.use, Animal_ID == levels(factor(determination.cohort$Animal_ID))[j])
    anim.spec.results <- table(anim.spec.data$WADDL.Movi.PCR)
    PropPos[j, i] <- anim.spec.results[1] / i
    TrueStatus[j] <- anim.spec.data$CarriageState[1]
  }
}

Det.Cohort.Tests <- as.data.frame(cbind(levels(factor(determination.cohort$Animal_ID)), PropPos, TrueStatus))
perc.cut <- .6

Det.Cohort.Tests$Correct1 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,2])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,2])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct2 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,3])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,3])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct3 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,4])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,4])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct4 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,5])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,5])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)

names(Det.Cohort.Tests) <- c("ID", "PropCorrect_PostTest1", "PropCorrect_PostTest2", "PropCorrect_PostTest3", "PropCorrect_PostTest4", "TrueStatus", "Correct1", "Correct2", "Correct3", "Correct4")
write.csv(Det.Cohort.Tests, "~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Documentation/BHSDetCohortTable.csv")

par(mfrow = c(1, 1), oma = c(2,0,0,0))
prop.correct <- c(sum(Det.Cohort.Tests$Correct1) / 12, sum(Det.Cohort.Tests$Correct2) / 12, sum(Det.Cohort.Tests$Correct3) / 12, sum(Det.Cohort.Tests$Correct4) / 12)
plot(prop.correct ~ c(1:4), xaxt = "n", pch = 16, ylim = c(0, 1), type = "b", ylab = "Accurate classification probability", xlab = "Number of tests")
axis(side = 1, at = c(1, 2, 3, 4), labels = c("1", "2", "3", "4"))
mtext(side = 1, line = 1, outer = T, cex = .9, adj = 0, "Classify an individual as chronic if proportion of positive tests is greater \n than .6")

#----------------------------------------------------------#
#-- bootstrap intervals for classification probabilities --#
#-- and numbers of tests ----------------------------------#
#----------------------------------------------------------#
nboot <- 1000
#-- approach is to resample all test results at individual level --#
test1.bootprops <- test1.bootprops <- test3.bootprops <- test4.bootprops <- rep(NA, nboot)
Resamp.Det.Cohort <- vector("list", nboot)
Resamp.prop.correct <- matrix(NA, nrow = nboot, ncol = 4)

for(i in 1:nboot){
  Resamp.Det.Cohort[[i]] <- Det.Cohort.Tests[sample(1:12, 12, rep = T), ]
  Resamp.prop.correct[i, ] <- c(sum(Resamp.Det.Cohort[[i]]$Correct1) / 12, sum(Resamp.Det.Cohort[[i]]$Correct2) / 12, sum(Resamp.Det.Cohort[[i]]$Correct3) / 12, sum(Resamp.Det.Cohort[[i]]$Correct4) / 12)
}

test1.lb <- quantile(Resamp.prop.correct[, 1], 0.025)
test1.ub <- quantile(Resamp.prop.correct[, 1], 0.975)
test2.lb <- quantile(Resamp.prop.correct[, 2], 0.025)
test2.ub <- quantile(Resamp.prop.correct[, 2], 0.975)
test3.lb <- quantile(Resamp.prop.correct[, 3], 0.025)
test3.ub <- quantile(Resamp.prop.correct[, 3], 0.975)
test4.lb <- quantile(Resamp.prop.correct[, 4], 0.025)
test4.ub <- quantile(Resamp.prop.correct[, 4], 0.975)
#-- END BOOTSTRAP --#

#--------------------------------------------------------------------------#
#-- Revised Accuracy post- number of tests plot with bootstrap intervals --#
#--------------------------------------------------------------------------#
require(graphics)
par(mfrow = c(1, 1), oma = c(1,0,0,0))
prop.correct <- c(sum(Det.Cohort.Tests$Correct1) / 12, sum(Det.Cohort.Tests$Correct2) / 12, sum(Det.Cohort.Tests$Correct3) / 12, sum(Det.Cohort.Tests$Correct4) / 12)
plot(prop.correct ~ c(1:4), xaxt = "n", pch = 16, ylim = c(0, 1), type = "b", ylab = "Accurate classification probability", xlab = "Number of tests")
segments(x0 = 1, x1 = 1, y0 = test1.lb, y1 = test1.ub, lty = 2, col = "grey40")
segments(x0 = 2, x1 = 2, y0 = test2.lb, y1 = test2.ub, lty = 2, col = "grey40")
segments(x0 = 3, x1 = 3, y0 = test3.lb, y1 = test3.ub, lty = 2, col = "grey40")
segments(x0 = 4, x1 = 4, y0 = test4.lb, y1 = test4.ub, lty = 2, col = "grey40")
axis(side = 1, at = c(1, 2, 3, 4), labels = c("1", "2", "3", "4"))

mtext(side = 1, line = 3, outer = T, cex = .9, adj = 0, "Classify an individual as chronic if proportion of positive tests is greater \nthan .6. Dashed grey lines are 95% bootstrapped confidence intervals for \nthe proportion of positive tests. No lines appear for tests 3 and 4, since all \nanimals in the determination cohort were correctly classified by that point.")

#-- END REVISED Accuracy post-number of tests plot --#

#---------------------#
#-- ELISA questions --#
#---------------------#

#-- 1) How does ELISA scale with PCR status for BHS?
#-- 2) Is ELISA inhibition level consistent within an animal through time?
#-- 3) Is inhibition level consistent with or predictive of infection/shedding state?

ManyAnimalSamples.full <- subset(data, Animal_ID %in% names(Animal.tab[which(Animal.tab >= 4)]))
ManyAnimalSamples <- subset(ManyAnimalSamples.full, Movipneumonia.ELISA != 0)

par(cex.axis = 1, las = 2, adj = .5)
plot(ManyAnimalSamples$Movipneumonia.ELISA ~ factor(ManyAnimalSamples$Animal_ID), type = "p", col = c("grey80", "grey80", rep("grey20", 7), "white", "grey20", "grey80"), ylab = "Movi ELISA", xlab = "")
leg.text <- c("not chronic", "chronic")
legend("bottomleft", bty = "n", leg.text, fill = c("grey20", "grey80"))
mtext(side = 1, line = 4, las = 1, "Animal ID")

#-- ELISA inhibition by PCR on a per-test (as opposed to within-sheep) basis --#
data.nozeros <- subset(data, Movipneumonia.ELISA != 0)
data.detindet <- subset(data.nozeros, WADDL.Movi.PCR != "Not detected")
quantile(data.detindet$Movipneumonia.ELISA, 0.025)
#-- 2.5th quantile is at ELISA = 49
data.notdetected <- subset(data.nozeros, WADDL.Movi.PCR == "Not detected")
table(data.notdetected$Movipneumonia.ELISA <= 49)
#-- 25/103 = 24% not detecteds are below the 2.5th ELISA quantile for detected or indeterminates across, Lostine, Asotin, and Black Butte

lostine <- subset(data.nozeros, Herd == "Lostine")
lostine.detindet <- subset(lostine, WADDL.Movi.PCR != "Not detected")
quantile(lostine.detindet$Movipneumonia.ELISA, 0.025)
#-- 2.5th quantile is at ELISA = 47.689

lostine.notdetected <- subset(lostine, WADDL.Movi.PCR == "Not detected")
table(lostine.notdetected$Movipneumonia.ELISA <= 47.689)
#-- 12/50 = 24% below the 2.5th ELISA quantile for detected or indeterminates at the Lostine

par(oma = c(2, 0, 0, 0))
plot(lostine$Movipneumonia.ELISA ~ jitter(as.numeric(lostine$WADDL.Movi.PCR), .25), xaxt = "n", xlab = "WADDL Movi PCR status", ylab = "Movi ELISA", ylim = c(0, 100))
abline(h = 49, col = "grey70", lty = 2)
axis(side = 1, at = c(1, 2, 3), labels = c("Detected", "Indeterminate", "Not dectected"))
mtext(side = 1, outer = T, line = 0, adj = 0, cex = .9, "All Lostine samples. Points are jittered on X-axis only (Movi ELISA values are \nexact).")

par(oma = c(2, 0, 0, 0))
plot(data.nozeros$Movipneumonia.ELISA ~ jitter(as.numeric(data.nozeros$WADDL.Movi.PCR), .25), col = ifelse(data.nozeros$Herd == "Lostine", "black", ifelse(data.nozeros$Herd == "Asotin", "grey80", "red")), pch = ifelse(data.nozeros$Herd == "Lostine", 1, 16), xaxt = "n", xlab = "WADDL Movi PCR status", ylab = "Movi ELISA", ylim = c(0, 100))
abline(h = 49, col = "grey70", lty = 2)
leg.text <- c("Lostine", "Asotin", "Black Butte")
legend("bottomleft", bty = "n", leg.text, pch = c(1, 16, 16), col = c("black", "grey80", "red"))
axis(side = 1, at = c(1, 2, 3), labels = c("Detected", "Indeterminate", "Not dectected"))
mtext(side = 1, outer = T, line = 0, adj = 0, cex = .9, "Lostine, Asotin, and Black Butte samples. Points are jittered on X-axis only (Movi \nELISA values are exact).")

#------------------------------------------------------------#
#-- Prediction of 4th test result for animals with 3 tests --#
#------------------------------------------------------------#
ThreeTests <- subset(data, Animal_ID %in% names(Animal.tab[which(Animal.tab == 3)]))
ThreeTestPosProp <- PredictedStatus <- rep(NA, length(names(Animal.tab[which(Animal.tab == 3)])))
for(i in 1:length(ThreeTestPosProp)){
  k <- subset(data, Animal_ID == names(Animal.tab[which(Animal.tab == 3)])[i])
  ThreeTestPosProp[i] <- ifelse(is.na(table(k$WADDL.Movi.PCR == "Detected")["TRUE"]) == T, 0, table(k$WADDL.Movi.PCR == "Detected")["TRUE"])
  PredictedStatus[i] <- ifelse(ThreeTestPosProp[i] >= 2, "chronic", "not chronic")
}

ThreeTestDataframe <- as.data.frame(cbind(names(Animal.tab[which(Animal.tab == 3)]), ThreeTestPosProp, PredictedStatus))
names(ThreeTestDataframe) <- c("ID", "Postive Tests", "Predicted Status")
write.csv(ThreeTestDataframe, "~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Documentation/BHSThreeTestStatusPrediction_18July2014.csv")

#----------------------------------------------------------#
#-- Prediction for final status for animals with 2 tests --#
#----------------------------------------------------------#

#-- question: what is the probability that an animal with 1 pos and 1 neg test winds up being chronic?  --#
TwoTests <- subset(data, Animal_ID %in% names(Animal.tab[which(Animal.tab == 2)]))
TwoTestTab <- table(factor(TwoTests$Animal_ID), TwoTests$WADDL.Movi.PCR)
TwoTestPred <- rep(NA, dim(TwoTestTab)[1])
for(i in 1:dim(TwoTestTab)[1]){
  TwoTestPred[i] <- ifelse(TwoTestTab[i, 1] == 2, "Chronic", ifelse(TwoTestTab[i, 3] == 2, "Not Chronic", "Indeterminate"))
}

TwoTestData <- cbind(TwoTestTab, TwoTestPred)

write.csv(TwoTestData, "~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Documentation/BHSTwoTestPredications_18July2014.csv")

#---------------------------------------------------#
#-- Tom's PCR data from 11 July 2014 ---------------#
#---------------------------------------------------#
#require(ROCR)
pcr.data <- read.csv("~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Data/BesserData_11July2014.csv", header = T)

#-- goal: model consistent positives (Result1, Result2) ~ ct_direct

#-- 1) Build response variable. use 2 positives as classification criterion for "consistent" positives. --#
pcr.data$May13Feb14 <- pcr.data$May13Feb14.binary <-rep(NA, dim(pcr.data)[1])
for(i in 1:dim(pcr.data)[1]){
  pcr.data$May13Feb14[i] <- ifelse(pcr.data$Result1[i] == "-" | pcr.data$Result2[i] == "-", NA, ifelse(pcr.data$Result1[i] == 1 & pcr.data$Result2[i] == 1, 1, 0))
}
pcr.data$May13Feb14.binary <- ifelse(pcr.data$May13Feb14 == 2 | pcr.data$May13Feb14 == 1, 1, 0)

pcr.data$consistent.pos <- rep(NA, dim(pcr.data)[1])
for(i in 1:dim(pcr.data)[1]){
  pcr.data$consistent.pos[i] <- pcr.data$May13Feb14[i]
}

#-- 2) track animals that were eliminated from analysis --#
complete.pcr.data.full <- subset(pcr.data, is.na(consistent.pos) == F)
  #-- 106 sheep sampled in both event 1 and event 2
table(complete.pcr.data.full$Result1, complete.pcr.data.full$Result2)
complete.pcr.data.noindet <- subset(complete.pcr.data.full, Result1 != 2 & Result2 != 2)
  #-- 92 sheep with no 2's
table(complete.pcr.data.noindet$Result1, complete.pcr.data.noindet$Result2)
table(is.na(complete.pcr.data.noindet$Ct_Direct))
  #-- 85 of those 92 also have Ct-Direct measures for event 3.
complete.pcr.data <- subset(complete.pcr.data.noindet, is.na(Ct_Direct) == F)

#-- 3) fit logistic regression model (with 0/0 animals included) --#
direct.model <- glm(complete.pcr.data$consistent.pos ~ complete.pcr.data$Ct_Direct, family = "binomial")

#-- 4) follow up with ROC analysis--#
complete.pcr.data$direct.preds <- predict(direct.model, complete.pcr.data, type = "response")
direct.labels <- complete.pcr.data$consistent.pos
#table(complete.pcr.data$direct.preds >= .5, complete.pcr.data$consistent.pos)
rocr.prediction <- prediction(complete.pcr.data$direct.preds, direct.labels)
rocr.performance <- performance(rocr.prediction, measure = "tpr", x.measure = "fpr")
auc.performance <- performance(rocr.prediction, measure = "auc")
  #-- AUC = 0.697 for direct as predictor of Result 1 & 2
#plot(rocr.performance)

#-- 5) Model that eliminates 2 0/0 results --#
complete.no.negs <- subset(complete.pcr.data, !(Result1 == 0 & Result2 == 0))

direct.model.nonegs <- glm(complete.no.negs$consistent.pos ~ complete.no.negs$Ct_Direct, family = "binomial")
complete.no.negs$direct.nonegs.preds <- predict(direct.model.nonegs, complete.no.negs, type = "response")
rocr.prediction.nonegs <- prediction(complete.no.negs$direct.nonegs.preds, complete.no.negs$consistent.pos)
rocr.performance.nonegs <- performance(rocr.prediction.nonegs, measure = "tpr", x.measure = "fpr")
auc.performance.nonegs <- performance(rocr.prediction.nonegs, measure = "auc")
#-- AUC is 0.702
#plot(rocr.performance.nonegs)
neg.negs <- subset(complete.pcr.data, Result1 == 0 & Result2 == 0)

#-- 6) Model 3: No Ct Direct in excess of 39 --#
complete.geq39 <- subset(complete.no.negs, Ct_Direct >= 39)
complete.sub39 <- subset(complete.no.negs, Ct_Direct <= 38.5)
direct.model.sub39 <- glm(complete.sub39$consistent.pos ~ complete.sub39$Ct_Direct, family = "binomial")
complete.sub39$direct.sub39.preds <- predict(direct.model.sub39, complete.sub39, type = "response")
rocr.prediction.sub39 <- prediction(complete.sub39$direct.sub39.preds, complete.sub39$consistent.pos)
rocr.performance.sub39 <- performance(rocr.prediction.sub39, measure = "tpr", x.measure = "fpr")
auc.performance.sub39 <- performance(rocr.prediction.sub39, measure = "auc")
#-- AUC is 0.744
plot(rocr.performance.nonegs)

#-- 7) plot of predicted probabilities of consistent carriage under each Model --#
complete.pcr.data$direct.preds <- predict(direct.model, complete.pcr.data, type = "response")
new.preds <- predict(direct.model, complete.pcr.data, type = "response", se.fit = T)
complete.pcr.data$se.preds <- new.preds$se.fit

complete.no.negs$direct.preds <- predict(direct.model.nonegs, complete.no.negs, type = "response")
new.preds.nonegs <- predict(direct.model.nonegs, complete.no.negs, type = "response", se.fit = T)
complete.no.negs$se.preds <- new.preds.nonegs$se.fit

complete.sub39$direct.preds <- predict(direct.model.sub39, complete.sub39, type = "response")
new.preds <- predict(direct.model.sub39, complete.sub39, type = "response", se.fit = T)
complete.sub39$se.preds <- new.preds$se.fit

par(mfrow = c(1, 1), oma = c(0, 1, 0, 0))
plot(complete.pcr.data$direct.preds ~ complete.pcr.data$Ct_Direct, type = "p", xlab = "Ct Direct", ylab = "Predicted probability of being consistently positive")
points(complete.pcr.data$direct.preds + (1.95 * complete.pcr.data$se.preds) ~ complete.pcr.data$Ct_Direct, col = "red", cex = .8, pch = 16)
points(complete.pcr.data$direct.preds - (1.95 * complete.pcr.data$se.preds) ~ complete.pcr.data$Ct_Direct, col = "red", cex = .8, pch = 16)

plot(complete.no.negs$direct.preds ~ complete.no.negs$Ct_Direct, type = "p", xlab = "Ct Direct", ylab = "Predicted probability of being consistently positive")
points(complete.no.negs$direct.preds + (1.95 * complete.no.negs$se.preds) ~ complete.no.negs$Ct_Direct, col = "red", cex = .8, pch = 16)
points(complete.no.negs$direct.preds - (1.95 * complete.no.negs$se.preds) ~ complete.no.negs$Ct_Direct, col = "red", cex = .8, pch = 16)

plot(complete.sub39$direct.preds ~ complete.sub39$Ct_Direct, type = "p", xlab = "Ct Direct", ylab = "Predicted probability of being consistently positive")
points(complete.sub39$direct.preds + (1.95 * complete.sub39$se.preds) ~ complete.sub39$Ct_Direct, col = "red", cex = .8, pch = 16)
points(complete.sub39$direct.preds - (1.95 * complete.sub39$se.preds) ~ complete.sub39$Ct_Direct, col = "red", cex = .8, pch = 16)

mtext(side = 1, outer = T, line = 8, adj = 0, cex = .9, "Probability of being consistently positive as estimated by a logistic regression model that contrasts individuals with two \n positive results in sampling events 1 and 2 to the same individual's Ct_Direct result from sampling event 3. Panel 1 shows fits \n with all 98 animals who had results for events 1 and 2, and Ct_Direct results for 3. Panel 2 shows the same model, fit after \n eliminating two individuals who were consistently negative. Panel 3 shows fits from the same model, but using a dataset \n that excludes both the negative-negative animals and 7 animals with Ct_Direct scores in excess of 38.5. Black dots are the \n point-estimates of probabilities of consistent shedding corresponding to each Ct_Direct value; red lines are 95% confidence \nlimits on those estimates.")



#---------------------------------------------------#
#-- composite model with both qPCR and ELISA -------#
#---------------------------------------------------#
# in preliminary models for BHS, qPCR is binary; needs to be continuous in future models. --#
# use data with known states: ManyAnimalSamples.full

# 1) label animals with a "true" state based on at least four tests. 
Chronic <- c("01LO34", "04LO12", "99L09")
ManyAnimalSamples.full$CarriageState <- ifelse(ManyAnimalSamples.full$Animal_ID %in% Chronic, 1, 0)

composite.bhs.glm <- glm(as.numeric(CarriageState) ~ Movipneumonia.ELISA + WADDL.Movi.PCR, data = ManyAnimalSamples.full, family = "binomial")
composite.bhs.glm.inter <- glm(as.numeric(CarriageState) ~ Movipneumonia.ELISA * WADDL.Movi.PCR, data = ManyAnimalSamples.full, family = "binomial")
anova(composite.bhs.glm, composite.bhs.glm.inter)

summary(composite.bhs.glm)
summary(composite.bhs.glm.inter)
#-- incorporate REs --#
require(lme4)
composite.bhs.glmer <- glmer(as.numeric(CarriageState) ~ Movipneumonia.ELISA + WADDL.Movi.PCR + (1 | Animal_ID) + (1|Herd), data = ManyAnimalSamples.full, family = "binomial")

# clustering experiment #
require(cluster)
daisy.data <- subset(data, select = c("Movipneumonia.ELISA", "WADDL.Movi.PCR"))
class(daisy.data$WADDL.Movi.PCR)
class(daisy.data$Movipneumonia.ELISA)

daisy.test <- daisy(daisy.data)
agnes.test <- agnes(daisy.test)
agnes.leafnames.ordered <- data$Animal_ID[as.numeric(as.character(agnes.test$order.lab))]
agnes.leafcolors <- ifelse(as.character(agnes.leafnames.ordered) %in% levels(factor(ManyAnimalSamples.full$Animal_ID)), ifelse(as.character(agnes.leafnames.ordered) %in% Chronic, "red", "blue"), "white")
agnes.dendro <- as.dendrogram(agnes.test)
plot(agnes.dendro,  nodePar = list(lab.col = agnes.leafcolors))

# install.packages("dendextend")
require(dendextend)
labels_colors(agnes.dendro) <- agnes.leafcolors
labels <- as.character(agnes.leafnames.ordered)
plot(agnes.dendro)

# mixture discriminant analysis experiment #
