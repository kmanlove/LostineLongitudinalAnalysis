full.data <- read.csv("./Data/LostineSample2016Final_WithClasses_20160927.csv", header = T)
full.data <- subset(full.data, cap_bioyr >= 2008)
full.data <- subset(full.data, movi_qpcr %in% c("Detected", "Not detected"))
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")

# remove all animals with only one observation
observed.once <- names(table(full.data$id))[which(table(full.data$id) == 1)]

data <- subset(full.data, !(id %in% observed.once) & age_class == "Adult")
data$id <- factor(data$id)
data$qPCRResult <- ifelse(data$movi_qpcr == "Detected" | data$movi_qpcr == "Indeterminate", "P", "N")

individ.bioyr <- unique(subset(data, select = c(id, cap_bioyr)))

# build list of all tests for each animal within a biological year
individ.bioyr.list <- vector("list", dim(individ.bioyr)[1])
number.tests.per.individ.bioyr <- rep(NA, dim(individ.bioyr)[1])

for(i in 1:length(individ.bioyr.list)){
  individ.bioyr.list[[i]] <- subset(data, id == individ.bioyr[i, 1] & cap_bioyr == individ.bioyr[i, 2])
  number.tests.per.individ.bioyr[i] <- dim(individ.bioyr.list[[i]])[1]
}

individ.bioyr.list.reduced <- individ.bioyr.list[number.tests.per.individ.bioyr != 1]
number.tests <- tests.pos <- tests.neg <- genotype.vsn <- rep(NA, length(individ.bioyr.list.reduced))
for(i in 1:length(individ.bioyr.list.reduced)){
  number.tests[i] <- dim(individ.bioyr.list.reduced[[i]])[1]
  tests.pos[i] <- dim(subset(individ.bioyr.list.reduced[[i]], movi_qpcr == "Detected"))[1]
  tests.neg[i] <- dim(subset(individ.bioyr.list.reduced[[i]], movi_qpcr == "Not detected"))[1]
  genotype.vsn[i] <- paste(individ.bioyr.list.reduced[[i]]$qPCRResult[1], "/", individ.bioyr.list.reduced[[i]]$qPCRResult[2], sep = "")
}

# first test results
firsttests <- subset(full.data, FirstTestInd == 1 & age_class == "Adult")
table(firsttests$movi_qpcr)
prop.pos <- 18 / (18+30)
prop.neg <- 30/(18+30)

# Expected number with 2 pos if system at Hardy-Weinberg keq
# expected pos:
n.pos.pos <- (prop.pos ^ 2) * 41
n.neg.neg <- (prop.neg ^ 2) * 41
n.switch <- 2 * prop.pos * prop.neg * 41
expected <- c(n.neg.neg, n.switch, n.pos.pos)
observed <- table(tests.pos[number.tests == 2])
observed <- c(19, 9, 13)

obs.0.ci <- binom.test(x = observed[1], n = sum(observed))$conf.int * sum(observed)
obs.1.ci <- binom.test(observed[2], sum(observed))$conf.int * sum(observed)
obs.2.ci <- binom.test(observed[3], sum(observed))$conf.int * sum(observed)

exp.0.ci <- binom.test(x = ceiling(expected[1]), n = sum(observed))$conf.int * sum(observed)
exp.1.ci <- binom.test(ceiling(expected[2]), sum(observed))$conf.int * sum(observed)
exp.2.ci <- binom.test(ceiling(expected[3]), sum(observed))$conf.int * sum(observed)

require(hwde)
hw.test <- hwexact(obs.hom1 = 19, obs.hets = 9, obs.hom2 = 13)

#---------------------------#
#-- OLD VSN ----------------#
#---------------------------#
# sum all test results to get overall distribution
prop.pos <- sum(tests.pos) / sum(number.tests)
prop.neg <- sum(tests.neg) / sum(number.tests)
sum(table(tests.pos[number.tests == 2]))
table(tests.pos[number.tests == 3])

# Expected number with 2 pos if system at Hardy-Weinberg keq
# expected pos:
n.pos.pos <- (prop.pos ^ 2) * 36
n.neg.neg <- (prop.neg ^ 2) * 36
n.switch <- 2 * prop.pos * prop.neg * 36
expected <- c(n.neg.neg, n.switch, n.pos.pos)
observed <- table(tests.pos[number.tests == 2])

obs.0.ci <- binom.test(x = observed[1], n = sum(observed))$conf.int * sum(observed)
obs.1.ci <- binom.test(observed[2], sum(observed))$conf.int * sum(observed)
obs.2.ci <- binom.test(observed[3], sum(observed))$conf.int * sum(observed)

exp.0.ci <- binom.test(x = ceiling(expected[1]), n = sum(observed))$conf.int * sum(observed)
exp.1.ci <- binom.test(ceiling(expected[2]), sum(observed))$conf.int * sum(observed)
exp.2.ci <- binom.test(ceiling(expected[3]), sum(observed))$conf.int * sum(observed)


# build table for plotting
par(mfrow = c(1, 2))
hw.tab <- rbind(expected, observed)
barplot(t(hw.tab), beside = T, ylim = c(0, 30), las = 1, 
        ylab = "Frequency", xlab = "Test Results",
        names.arg = c("Expected under H-W", "Observed"),
        legend.text = c("Never positive", "Positive once", "Positive twice"),
        args.legend = list(y = 30, x = 4.5, bty = "n", cex = .8),
        main = "Within-year")
segments(x0 = 1.5, x1 = 1.5, y0 = exp.0.ci[1], y1 = exp.0.ci[2])
segments(x0 = 2.5, x1 = 2.5, y0 = exp.1.ci[1], y1 = exp.1.ci[2])
segments(x0 = 3.5, x1 = 3.5, y0 = exp.2.ci[1], y1 = exp.2.ci[2])

segments(x0 = 5.5, x1 = 5.5, y0 = obs.0.ci[1], y1 = obs.0.ci[2])
segments(x0 = 6.5, x1 = 6.5, y0 = obs.1.ci[1], y1 = obs.1.ci[2])
segments(x0 = 7.5, x1 = 7.5, y0 = obs.2.ci[1], y1 = obs.2.ci[2])


#---------------------------#
#-- X-year Hardy-Weinberg --#
#---------------------------#
# xyr <- subset(data, is.na(XyearStatus_2) == F)
# xyr.tab <- table(xyr$XyearStatus_2)

# prop.pos <- (23 * 2 + 12 + 11) / (2*(30 + 11 + 8 + 21))
# prop.neg <- (30 * 2 + 12 + 11) / (2*(30 + 11 + 8 + 21))
prop.pos.xyr <- (23 * 2 + 12 + 11) / (2*(30 + 11 + 8 + 21))
prop.neg.xyr <- (30 * 2 + 12 + 11) / (2*(30 + 11 + 8 + 21))

n.pos.pos.xyr <- (prop.pos.xyr ^ 2) * 76
n.neg.neg.xyr <- (prop.neg.xyr ^ 2) * 76
n.switch.xyr <- 2 * prop.pos.xyr * prop.neg.xyr * 76
expected.xyr <- c(n.neg.neg.xyr, n.switch.xyr, n.pos.pos.xyr)
observed.xyr <- c(30, 23, 23)

obs.0.ci.xyr <- binom.test(x = observed.xyr[1], n = sum(observed.xyr))$conf.int * sum(observed.xyr)
obs.1.ci.xyr <- binom.test(observed.xyr[2], sum(observed.xyr))$conf.int * sum(observed.xyr)
obs.2.ci.xyr <- binom.test(observed.xyr[3], sum(observed.xyr))$conf.int * sum(observed.xyr)

exp.0.ci.xyr <- binom.test(x = ceiling(expected.xyr[1]), n = sum(observed.xyr))$conf.int * sum(observed.xyr)
exp.1.ci.xyr <- binom.test(ceiling(expected.xyr[2]), sum(observed.xyr))$conf.int * sum(observed.xyr)
exp.2.ci.xyr <- binom.test(ceiling(expected.xyr[3]), sum(observed.xyr))$conf.int * sum(observed.xyr)

hw.tab.xyr <- rbind(expected.xyr, observed.xyr)
barplot(t(hw.tab.xyr), beside = T, ylim = c(0, 60), las = 1, 
        ylab = "Frequency", xlab = "Test Results",
        names.arg = c("Expected under H-W", "Observed"),
        legend.text = c("Never positive", "Positive once", "Positive twice"),
        args.legend = list(y = 60, x = 4.5, bty = "n", ncol = 1, cex = .8),
        main = "Between-year")
segments(x0 = 1.5, x1 = 1.5, y0 = exp.0.ci.xyr[1], y1 = exp.0.ci.xyr[2])
segments(x0 = 2.5, x1 = 2.5, y0 = exp.1.ci.xyr[1], y1 = exp.1.ci.xyr[2])
segments(x0 = 3.5, x1 = 3.5, y0 = exp.2.ci.xyr[1], y1 = exp.2.ci.xyr[2])

segments(x0 = 5.5, x1 = 5.5, y0 = obs.0.ci.xyr[1], y1 = obs.0.ci.xyr[2])
segments(x0 = 6.5, x1 = 6.5, y0 = obs.1.ci.xyr[1], y1 = obs.1.ci.xyr[2])
segments(x0 = 7.5, x1 = 7.5, y0 = obs.2.ci.xyr[1], y1 = obs.2.ci.xyr[2])


#-----------------------------#
#-- Full figure --------------#
#-----------------------------#
# svg("./Plots/HWPlots_20160623.svg", height = 5, width = 6)
par(mfrow = c(1, 2))
hw.tab <- rbind(expected, observed)
barplot(t(hw.tab), beside = T, ylim = c(0, 30), las = 1, 
        ylab = "Frequency", xlab = "Test Results",
        names.arg = c("Expected under H-W", "Observed"),
        legend.text = c("Never positive", "Positive once", "Positive twice"),
        args.legend = list(y = 30, x = 4.5, bty = "n", cex = .8),
        main = "Within-year")
segments(x0 = 1.5, x1 = 1.5, y0 = exp.0.ci[1], y1 = exp.0.ci[2])
segments(x0 = 2.5, x1 = 2.5, y0 = exp.1.ci[1], y1 = exp.1.ci[2])
segments(x0 = 3.5, x1 = 3.5, y0 = exp.2.ci[1], y1 = exp.2.ci[2])

segments(x0 = 5.5, x1 = 5.5, y0 = obs.0.ci[1], y1 = obs.0.ci[2])
segments(x0 = 6.5, x1 = 6.5, y0 = obs.1.ci[1], y1 = obs.1.ci[2])
segments(x0 = 7.5, x1 = 7.5, y0 = obs.2.ci[1], y1 = obs.2.ci[2])

barplot(t(hw.tab.xyr), beside = T, ylim = c(0, 60), las = 1, 
        ylab = "Frequency", xlab = "Test Results",
        names.arg = c("Expected under H-W", "Observed"),
        legend.text = c("Never positive", "Positive once", "Positive twice"),
        args.legend = list(y = 60, x = 4.5, bty = "n", ncol = 1, cex = .8),
        main = "Between-year")
segments(x0 = 1.5, x1 = 1.5, y0 = exp.0.ci.xyr[1], y1 = exp.0.ci.xyr[2])
segments(x0 = 2.5, x1 = 2.5, y0 = exp.1.ci.xyr[1], y1 = exp.1.ci.xyr[2])
segments(x0 = 3.5, x1 = 3.5, y0 = exp.2.ci.xyr[1], y1 = exp.2.ci.xyr[2])

segments(x0 = 5.5, x1 = 5.5, y0 = obs.0.ci.xyr[1], y1 = obs.0.ci.xyr[2])
segments(x0 = 6.5, x1 = 6.5, y0 = obs.1.ci.xyr[1], y1 = obs.1.ci.xyr[2])
segments(x0 = 7.5, x1 = 7.5, y0 = obs.2.ci.xyr[1], y1 = obs.2.ci.xyr[2])
# dev.off()