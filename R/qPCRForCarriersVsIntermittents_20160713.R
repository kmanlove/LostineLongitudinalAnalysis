samples <- read.csv("./Data/LostineSample2016Final.csv", header = T)

qpcrs <- subset(samples, is.na(movi_direct_pcr) == F)
qpcrs_testpos <- subset(samples, is.na(movi_direct_pcr) == F & movi_qpcr %in% c("Detected"))
qpcrs_testpos$carriage_twoPosOneYear <- factor(qpcrs_testpos$carriage_twoPosOneYear, 
                                               levels = c("carrier", "intermittent", "negative"))


par(mfrow = c(1, 2), las = 1, cex.axis = .8)
boxplot(qpcrs$movi_direct_pcr ~ as.numeric(as.factor(qpcrs$carriage_twoPosOneYear)), col = "grey60",
        xlim = c(0, 4), xaxt = "n", ylab = "qPCR cycle threshold",
        main = "All tests")
axis(side = 1, at = c(1, 2, 3), 
     labels = c("carrier", "intermittent", "negative"), cex = .6)
boxplot(qpcrs_testpos$movi_direct_pcr ~ as.numeric(qpcrs_testpos$carriage_twoPosOneYear), col = "grey60",
        xlim = c(0, 4), xaxt = "n", ylab = "qPCR cycle threshold",
        main = "Positive tests")
axis(side = 1, at = c(1, 2, 3), 
     labels = c("carrier", "intermittent", "negative"), cex = .6)


# test to compare qPCRs among test-pos carriers and test-pos intermittents
wilcox.testpos.out <- wilcox.test(qpcrs_testpos$movi_direct_pcr ~ qpcrs_testpos$carriage_twoPosOneYear)


# power analysis to examine differences
# effect might be two doublings shift...
alpha <- 0.05
# beta <- 0.80
# sample.size <- 
# variance <- 
# effect <- 
# effect.size <- effect / variance

# tabulate values on raw data
carrier.testpos <- subset(qpcrs_testpos, carriage_twoPosOneYear == "carrier")
inter.testpos <- subset(qpcrs_testpos, carriage_twoPosOneYear == "intermittent")
n.carrier <- dim(carrier.testpos)[1]
n.inter <- dim(inter.testpos)[1]
sd.carrier.ct <- sd(carrier.testpos$movi_direct_pcr)
sd.inter.ct <- sd(inter.testpos$movi_direct_pcr)
variance.in <- (sqrt(sd.carrier.ct/n.carrier + sd.inter.ct/n.inter))^2

# generate new data
baseline.ct <- 1
effects <- seq(0, 12, length.out = 25)
n.reps <- 2000
reject <- matrix(NA, nrow = length(effects), ncol = n.reps)
simmed.carriers <- simmed.inters <- wilcox.out <- vector("list", length(effects))
for(i in 1:length(effects)){
  simmed.carriers[[i]] <- simmed.inters[[i]] <- vector("list", n.reps)
  for(j in 1:n.reps){
    simmed.carriers[[i]][[j]] <- rnorm(mean = baseline.ct, sd = sd.carrier.ct, n = n.carrier)
    simmed.inters[[i]][[j]] <- rnorm(mean = baseline.ct - effects[i], sd = sd.inter.ct, n = n.inter)
    simmed.cts <- c(simmed.carriers[[i]][[j]], simmed.inters[[i]][[j]])
    simmed.states <- c(rep("carrier", n.carrier), rep("inter", n.inter))
    wilcox.out[[i]][[j]] <- wilcox.test(simmed.cts ~ simmed.states)
    reject[i, j] <- ifelse(wilcox.out[[i]][[j]]$p.value <= alpha, 1, 0)
  }
  print(i)
}

# table rejection rates for each effect size
rejection.prop <- matrix(NA, ncol = 2, nrow = length(effects))
for(i in 1:length(effects)){
  rejection.prop[i, ] <- table(reject[i, ])/dim(reject)[2]
}

# "1"s are correct rejections
par(mfrow = c(1, 1))
plot(rejection.prop[,2] ~ effects, 
     ylab = "P(Reject null | Alternative true)", 
     xlab = "True difference in cts", type = "l")
abline(h = .8, lty = 2, col = "red")
abline(v = 9.5, lty = 2, col = "red")
