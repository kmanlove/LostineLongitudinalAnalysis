full.data <- read.csv("./Data/LostineSample2016Final_WithClasses_20160927.csv", header = T)
full.data <- subset(full.data, movi_qpcr %in% c("Detected", "Not detected"))
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")
full.data <- subset(full.data, cap_bioyr >= 2008)


# full.data <- read.csv("./Data/LostineSample2016Final.csv", header = T)
# full.data$id <- full.data$Animal.ID
# full.data$movi_qpcr <- full.data$WADDL.Movi.PCR
# full.data$cap_bioyr <- full.data$CAP_BIOYR
# full.data <- subset(full.data, movi_qpcr != "Indeterminate")
# full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")

full.data <- subset(full.data, 
                    Age.Class != "Lamb" & Age.Class != "Yearling" &
                    !(movi_qpcr %in% c("", "Indeterminate")))

table(factor(full.data$id))
# animals with 3 or more samples
many.samples <- names(table(factor(full.data$id)))[which(table(factor(full.data$id)) >= 3)]

full.data <- subset(full.data, id %in% many.samples)
full.data$id <- factor(full.data$id)
full.data$movi_qpcr <- factor(full.data$movi_qpcr)

# table by test outcome
movi.tab <- table(full.data$id, full.data$movi_qpcr)
prop.pos <- movi.tab[,1] / (movi.tab[,1] + movi.tab[,2])

#chronics <- names(prop.pos)[which(prop.pos >= 0.787)]
chronics <- levels(factor(subset(full.data, carrier_posConsecYears == "carrier")$id))

# bring in movi survival data
#movi.surv <- read.csv("./Data/MoviSurvivalDate_20160725.csv", header = T)
movi.surv <- read.csv("./Data/MoviSurvivalData_20160822.csv", header = T)
movi.surv$stopdate.Date <- as.Date(movi.surv$stopdate, format = "%m/%d/%Y")
movi.surv$startdate.Date <- as.Date(movi.surv$startdate, format = "%m/%d/%Y")

movi.surv$duration <- rep(NA, dim(movi.surv)[1])
for(i in 1:dim(movi.surv)[1]){
  movi.surv$duration[i] <- difftime(movi.surv$stopdate.Date[i], 
                                    movi.surv$startdate.Date[i])
}

movi.surv$CensoringStatus <- rep(NA, dim(movi.surv)[1])
for(i in 1:dim(movi.surv)[1]){
  if(movi.surv$rightcensor[i] == 1){
    if(movi.surv$leftcensor[i] == 1){
      movi.surv$CensoringStatus[i] <- "rightleft"
    } else if(movi.surv$leftcensor[i] == 0){
      movi.surv$CensoringStatus[i] <- "right"
    }
  } else if(movi.surv$leftcensor[i] == 1){
    movi.surv$CensoringStatus[i] <- "left"
  } else movi.surv$CensoringStatus[i] <- "none"
}

boxplot(movi.surv$duration ~ movi.surv$CensoringStatus,
        col = "grey80", ylim = c(0, 1200), 
        las = 1, xlab = "Censoring status", ylab = "Duration")
text(y = rep(1200, 4), x = c(1:4), cex = .8, 
     labels = c("N = 16", "N = 6", "N = 13", "N = 15"))

# extract infections in chronic animals
chronic.movi.surv <- subset(movi.surv, ID %in% chronics)
hist(movi.surv$duration, col = "grey80", breaks = 20)
dim(chronic.movi.surv) # 11 infections in "chronic" animals
length(levels(factor(chronic.movi.surv$ID))) # 10 animals are "chronic"

mean(chronic.movi.surv$duration)
median(chronic.movi.surv$duration)
quantile(chronic.movi.surv$duration, c(.05, .25, .5, .75, .95))

# long infections in "transient" animals
transient.movi.surv <- subset(movi.surv, !(ID %in% chronics))
long.transient.rows <- which(transient.movi.surv$duration >= median(chronic.movi.surv$duration))
short.transient.rows <- which(transient.movi.surv$duration < median(chronic.movi.surv$duration))

long.transients <- transient.movi.surv[long.transient.rows, ]
short.transients <- transient.movi.surv[short.transient.rows, ]

table(long.transients$rightcensor)
table(short.transients$rightcensor)
table(long.transients$rightcensor == 0, long.transients$leftcensor == 0)
table(short.transients$rightcensor == 0, short.transients$leftcensor == 0)
table(chronic.movi.surv$rightcensor == 0, chronic.movi.surv$leftcensor == 0)


par(mfrow = c(1, 2))
boxplot(transient.movi.surv$duration ~ transient.movi.surv$rightcensor, col = "grey80",
        names = c("not \ncensored", "censored"), 
        las = 1,  ylab = "infection duration",
        ylim = c(0, 1200),
        main = "infection duration in \ntransient animals")
abline(h = 400, lwd = .5, lty = 2)
abline(h = 800, lwd = .5, lty = 2)
boxplot(chronic.movi.surv$duration ~ chronic.movi.surv$rightcensor, col = "grey80",
        names = c("not \ncensored", "censored"), 
        las = 1,  ylab = "infection duration",
        ylim = c(0, 1200),
        main = "infection duration in \nchronic animals")
abline(h = 400, lwd = .5, lty = 2)
abline(h = 800, lwd = .5, lty = 2)


par(mfrow = c(1, 3))
boxplot(chronic.movi.surv$duration ~ chronic.movi.surv$rightcensor, col = "grey80",
        names = c("not censored", "censored"), 
        las = 1,  ylab = "infection duration",
        ylim = c(0, 1200),
        main = "infection duration in \nchronic animals")
boxplot(long.transients$duration ~ long.transients$rightcensor, col = "grey80",
        names = c("not censored", "censored"), las = 1,  ylab = "infection duration",
        main = "long infections in \ntransient animals",
        ylim = c(0, 1200))
boxplot(short.transients$duration ~ short.transients$rightcensor, col = "grey80",
        names = c("not censored", "censored"), las = 1,  ylab = "infection duration",
        main = "short infections in \ntransient animals",
        ylim = c(0, 1200))

par(mfrow = c(1, 3))
plot(movi.surv$duration ~ as.numeric(factor(movi.surv$CensoringStatus, levels = levels(factor(movi.surv$CensoringStatus)))),
        col = rgb(.3, .3, .3, alpha = .6), pch = 16, ylim = c(0, 1200), 
        las = 1, xlab = "Censoring status", main = "All infections",
     ylab = "Duration", cex = 2, xaxt = "n", xlim = c(.5, 4.5))
text(y = rep(1200, 4), x = c(1:4), cex = .8, 
     labels = c("N = 16", "N = 6", "N = 13", "N = 15"))
axis(side = 1, at = c(1:4), labels = levels(factor(movi.surv$CensoringStatus)))
plot(chronic.movi.surv$duration ~ as.numeric(factor(chronic.movi.surv$CensoringStatus, levels = levels(factor(movi.surv$CensoringStatus)))), 
     col = rgb(.3, .3, .3, alpha = .6), pch = 16, ylim = c(0, 1200), 
     las = 1,  ylab = "infection duration",
     main = "infection duration in \nchronic animals (pos >78% of tests)",
     cex = 2, xlab = "Censoring status", 
     xaxt = "n", xlim = c(.5, 4.5))
text(y = rep(1200, 4), x = c(1:4), cex = .8, 
     labels = c("N = 1", "N = 0", "N = 1", "N = 9"))
axis(side = 1, at = c(1:4), labels = levels(factor(movi.surv$CensoringStatus)))
plot(transient.movi.surv$duration ~ as.numeric(factor(transient.movi.surv$CensoringStatus, levels = levels(factor(movi.surv$CensoringStatus)))), 
     las = 1,  ylab = "infection duration",
     col = rgb(.3, .3, .3, alpha = .6), pch = 16, 
     ylim = c(0, 1200), xlab = "Censoring status",
     main = "infection duration in \ntransient animals", 
     cex = 2, xaxt = "n", xlim = c(.5, 4.5))
text(y = rep(1200, 4), x = c(1:4), cex = .8, 
     labels = c("N = 15", "N = 6", "N = 12", "N = 6"))
axis(side = 1, at = c(1:4), labels = levels(factor(movi.surv$CensoringStatus)))
