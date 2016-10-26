#--------------------------------------------------------#
#-- Frances' longitudinal data from the Lostine/Asotin --#
#--------------------------------------------------------#
require(graphics)

data <- read.csv("./Data/Old/LongSamplingData_05July2014.csv", header = T, sep = "\t")
data <- read.csv("./Data/LostineMoviPastLung_150630.csv", header = T)
full.data <- read.csv("./Data/LostineSample2016Final.csv", header = T)
head(data)

# data$WADDL.Movi.PCR <- ifelse(data$WADDL.Movi.PCR == "Indeterminant", "Indeterminate", as.character(data$WADDL.Movi.PCR))
# data$WADDL.Movi.PCR <- factor(data$WADDL.Movi.PCR)
# table(data$WADDL.Movi.PCR)

#-- build and explore determination cohort --#
# Animal.tab <- table(data$Animal_ID)
# AnimalByOutcome.tab <- table(data$Animal_ID, data$WADDL.Movi.PCR)
# Animal.tab[which(Animal.tab >= 4)]
# Animal.tab[which(Animal.tab == 3)]

Animal.tab <- table(data$id)
AnimalByOutcome.tab <- table(data$id, data$movi_qpcr)
Animal.tab[which(Animal.tab >= 4)]
Animal.tab[which(Animal.tab == 3)]
ageclass.tab <- table(data$id, data$age_class)
ad.ids <- row.names(ageclass.tab[which(ageclass.tab[, 2] == 0 & ageclass.tab[ ,3] == 0), ])

ManySamplesOutcomes.tab <- AnimalByOutcome.tab[which(Animal.tab >= 4), ]
ManySamplesProps <- ManySamplesOutcomes.tab[ ,1] / (ManySamplesOutcomes.tab[ ,1] + ManySamplesOutcomes.tab[ ,3])

# determination.cohort <- subset(data, Animal_ID %in% names(which(Animal.tab >= 4)))
# Chronic <- c("01LO34", "04LO12", "99L09")
# determination.cohort$CarriageState <- ifelse(determination.cohort$Animal_ID %in% Chronic, "chronic", "not.chronic")
determination.cohort <- subset(data, id %in% names(which(Animal.tab >= 4)) & id %in% ad.ids)
status.tab <- table(determination.cohort$movi_qpcr, factor(determination.cohort$id))
# to be chronic, need 75% of test results positive. 
Chronic <- c("01LO32", "01LO34", "02LO43", "03LO48", "04LO12", "04LO74", "99L09")
determination.cohort$CarriageState <- ifelse(determination.cohort$id %in% Chronic, "chronic", "not.chronic")

# PropPos <- matrix(NA, ncol = 4, nrow = length(levels(factor(determination.cohort$Animal_ID))))
# TrueStatus <- rep(NA, length(levels(factor(determination.cohort$Animal_ID))))
PropPos <- matrix(NA, ncol = 7, nrow = length(levels(factor(determination.cohort$id))))
TrueStatus <- rep(NA, length(levels(factor(determination.cohort$id))))
# for(i in 1:4){
#   tests.to.use <- subset(determination.cohort, TestNumber <= i)
#   for(j in 1:length(levels(factor(determination.cohort$Animal_ID)))){
#     anim.spec.data <- subset(tests.to.use, Animal_ID == levels(factor(determination.cohort$Animal_ID))[j])
#     anim.spec.results <- table(anim.spec.data$WADDL.Movi.PCR)
#     PropPos[j, i] <- anim.spec.results[1] / i
#     TrueStatus[j] <- anim.spec.data$CarriageState[1]
#   }
# }
for(i in 1:7){
  tests.to.use <- subset(determination.cohort, TestNumber <= i)
  for(j in 1:length(levels(factor(determination.cohort$id)))){
    anim.spec.data <- subset(tests.to.use, id == levels(factor(determination.cohort$id))[j])
    anim.spec.results <- table(anim.spec.data$movi_qpcr)
    PropPos[j, i] <- anim.spec.results[1] / i
    TrueStatus[j] <- anim.spec.data$CarriageState[1]
  }
}

Det.Cohort.Tests <- as.data.frame(cbind(levels(factor(determination.cohort$id)), PropPos, TrueStatus))
perc.cut <- .6

Det.Cohort.Tests$Correct1 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,2])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,2])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct2 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,3])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,3])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct3 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,4])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,4])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct4 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,5])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,5])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct5 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,6])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,6])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct6 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,7])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,7])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)
Det.Cohort.Tests$Correct7 <- ifelse(as.numeric(as.character(Det.Cohort.Tests[ ,8])) >= perc.cut & Det.Cohort.Tests$TrueStatus == "chronic" | 
                                      as.numeric(as.character(Det.Cohort.Tests[ ,8])) <= perc.cut & Det.Cohort.Tests$TrueStatus == "not.chronic", 1, 0)

names(Det.Cohort.Tests) <- c("ID", "PropCorrect_PostTest1", "PropCorrect_PostTest2", 
                             "PropCorrect_PostTest3", 
                             "PropCorrect_PostTest4", 
                             "PropCorrect_PostTest5", 
                             "PropCorrect_PostTest6", 
                             "PropCorrect_PostTest7", 
                             "TrueStatus", 
                             "Correct1", 
                             "Correct2", 
                             "Correct3", 
                             "Correct4", 
                             "Correct5", 
                             "Correct6", 
                             "Correct7")
# write.csv(Det.Cohort.Tests, "~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Documentation/BHSDetCohortTable.csv")

#----------------------------------------------------------#
#-- bootstrap intervals for classification probabilities --#
#-- and numbers of tests ----------------------------------#
#----------------------------------------------------------#
nboot <- 1000
#-- approach is to resample all test results at individual level --#
test1.bootprops <- test1.bootprops <- test3.bootprops <- test4.bootprops <- test5.bootprops <- test6.bootprops <- test7.bootprops <- rep(NA, nboot)
Resamp.Det.Cohort <- vector("list", nboot)
Resamp.prop.correct <- matrix(NA, nrow = nboot, ncol = 7)
n.animals <- dim(status.tab)[2]

for(i in 1:nboot){
  Resamp.Det.Cohort[[i]] <- Det.Cohort.Tests[sample(1:n.animals, n.animals, rep = T), ]
  Resamp.prop.correct[i, ] <- c(sum(Resamp.Det.Cohort[[i]]$Correct1) / n.animals, 
                                sum(Resamp.Det.Cohort[[i]]$Correct2) / n.animals, 
                                sum(Resamp.Det.Cohort[[i]]$Correct3) / n.animals, 
                                sum(Resamp.Det.Cohort[[i]]$Correct4) / n.animals, 
                                sum(Resamp.Det.Cohort[[i]]$Correct5) / n.animals
#                                , 
#                                sum(Resamp.Det.Cohort[[i]]$Correct6) / n.animals, 
#                                sum(Resamp.Det.Cohort[[i]]$Correct7) / n.animals)
}

test1.lb <- quantile(Resamp.prop.correct[, 1], 0.025)
test1.ub <- quantile(Resamp.prop.correct[, 1], 0.975)
test2.lb <- quantile(Resamp.prop.correct[, 2], 0.025)
test2.ub <- quantile(Resamp.prop.correct[, 2], 0.975)
test3.lb <- quantile(Resamp.prop.correct[, 3], 0.025)
test3.ub <- quantile(Resamp.prop.correct[, 3], 0.975)
test4.lb <- quantile(Resamp.prop.correct[, 4], 0.025)
test4.ub <- quantile(Resamp.prop.correct[, 4], 0.975)
test5.lb <- quantile(Resamp.prop.correct[, 5], 0.025)
test5.ub <- quantile(Resamp.prop.correct[, 5], 0.975)
# test6.lb <- quantile(Resamp.prop.correct[, 6], 0.025)
# test6.ub <- quantile(Resamp.prop.correct[, 6], 0.975)
# test7.lb <- quantile(Resamp.prop.correct[, 7], 0.025)
# test7.ub <- quantile(Resamp.prop.correct[, 7], 0.975)
#-- END BOOTSTRAP --#

#--------------------------------------------------------------------------#
#-- Revised Accuracy post- number of tests plot with bootstrap intervals --#
#--------------------------------------------------------------------------#
par(mfrow = c(1, 1), oma = c(1,0,0,0))
prop.correct <- c(sum(Det.Cohort.Tests$Correct1) / n.animals, 
                  sum(Det.Cohort.Tests$Correct2) / n.animals, 
                  sum(Det.Cohort.Tests$Correct3) / n.animals, 
                  sum(Det.Cohort.Tests$Correct4) / n.animals, 
                  sum(Det.Cohort.Tests$Correct5) / n.animals
#                   , 
#                   sum(Det.Cohort.Tests$Correct6) / n.animals, 
#                   sum(Det.Cohort.Tests$Correct7) / n.animals
                  )
plot(prop.correct ~ c(1:5), xaxt = "n", pch = 16, ylim = c(0, 1), type = "b", 
     ylab = "Accurate classification probability", xlab = "Number of tests", las = 1)
segments(x0 = 1, x1 = 1, y0 = test1.lb, y1 = test1.ub, lty = 1, col = "grey30")
segments(x0 = 2, x1 = 2, y0 = test2.lb, y1 = test2.ub, lty = 1, col = "grey30")
segments(x0 = 3, x1 = 3, y0 = test3.lb, y1 = test3.ub, lty = 1, col = "grey30")
segments(x0 = 4, x1 = 4, y0 = test4.lb, y1 = test4.ub, lty = 1, col = "grey30")
segments(x0 = 5, x1 = 5, y0 = test5.lb, y1 = test5.ub, lty = 1, col = "grey30")
#segments(x0 = 6, x1 = 6, y0 = test6.lb, y1 = test6.ub, lty = 2, col = "grey40")
#segments(x0 = 7, x1 = 7, y0 = test7.lb, y1 = test7.ub, lty = 2, col = "grey40")
axis(side = 1, at = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5"))

mtext(side = 1, line = 3, outer = T, cex = .9, adj = 0, "Classify an individual as chronic if proportion of positive tests is greater \nthan .6. Dashed grey lines are 95% bootstrapped confidence intervals for \nthe proportion of positive tests. No lines appear for tests 3 and 4, since all \nanimals in the determination cohort were correctly classified by that point.")

#-- END REVISED Accuracy post-number of tests plot --#