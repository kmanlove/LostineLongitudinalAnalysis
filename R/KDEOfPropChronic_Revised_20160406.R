#full.data <- read.csv("./Data/LostineMoviPastLung_150630.csv", header = T)
#full.data <- read.csv("./Data/LostineSampleData_20160406_v2.csv", header = T)
#full.data <- read.csv("./Data/LostineSampleData_20160505_v1.csv", header = T)
full.data <- read.csv("./Data/LostineSample2016Final.csv", header = T)
full.data$id <- full.data$Animal.ID
full.data$movi_qpcr <- full.data$WADDL.Movi.PCR
full.data$cap_bioyr <- full.data$CAP_BIOYR
full.data <- subset(full.data, movi_qpcr != "Indeterminate")
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")

full.data <- subset(full.data, Age.Class != "Lamb" & Age.Class != "Yearling")

# remove all animals with only one observation
observed.once <- names(table(full.data$id))[which(table(full.data$id) == 1 )]
observed.twice <-  names(table(full.data$id))[which(table(full.data$id) == 2 )]
  
data <- subset(full.data)
               # , !(id %in% observed.once | id %in% observed.twice))
data$id <- factor(data$id)
data$qPCRResult <- ifelse(data$movi_qpcr == "Detected" | data$movi_qpcr == "Indeterminate", "P", "N")

pcr.by.id <- table(data$qPCRResult, data$id)
prop.pos <- prop.table(pcr.by.id, margin = 2)[1, ]
d <- density((1-prop.pos), 
             weights = as.vector(table(data$id))/sum(table(data$id)), 
             from = 0, to = 1,
             bw = .1)


# bootstrapped CI for minimum of d
nboot <- 10000
resamp.min <- rep(NA, nboot)
for(i in 1:nboot){
  resamp.dat <- pcr.by.id[ ,sample(1:31, rep = T)]
  resamp.prop.pos <- prop.table(resamp.dat, margin = 2)[1, ]
  resamp.d <- density((1-resamp.prop.pos), 
                      weights = as.vector(resamp.dat[1,] + resamp.dat[2, ])/sum(resamp.dat[1,] + resamp.dat[2, ]), 
                      from = 0, to = 1,
                      bw = .1)
  resamp.min[i] <- resamp.d$x[which.min(resamp.d$y)]
  print(i)
}

quantile(resamp.min, c(0.1, 0.5, 0.9))

# reduced to just animals with 3 or more samples
data <- subset(full.data , !(id %in% observed.once | id %in% observed.twice))
data$id <- factor(data$id)

pcr.by.id.full <- table(full.data$qPCRResult, full.data$id)
prop.pos.full <- prop.table(pcr.by.id.full, margin = 2)[1, ]

pcr.by.id <- table(data$qPCRResult, data$id)
prop.pos <- prop.table(pcr.by.id, margin = 2)[1, ]
rbind(pcr.by.id, 1-prop.pos)

#svg("./Plots/PropPositiveFrequencies_20160623.svg", height = 5, width = 6)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 4))
hist(1 - prop.pos, breaks = 10, col = rgb(0, 0, 0, alpha = .5),
     xlim = c(0, 1), ylim = c(0, 10),
     ylab = "# Animals", 
     xlab = "Proportion positive", 
     las = 1,
     main = "No lambs or yearlings")
leg.text <- "Animals with 3 or more tests"
legend(y = 9.5, x = 0, leg.text, fill = "grey50", bty = "n")

# par(new = T)
# hist(1 - prop.pos.full, breaks = 10, col = rgb(0, 0, 0, alpha = .2),
#      xlim = c(0, 1), ylim = c(0, 20),
#      main = "", ylab = "# Animals", 
#      xlab = "Proportion positive", 
#      las = 1)
# leg.text <- "Animals with 3 or more tests"
# legend(y = 9, x = 0, leg.text, fill = "grey50", bty = "n")

#             prop.pos)
par(new = T)
plot(d, xlim = c(0, 1), ylim = c(0, max(d$y)+.2), yaxt = "n", xlab = "", 
     main = "", ylab = "")
mtext(side = 4, line = 3, "")
axis(side = 4, at = seq(0, max(d$y), by = .25), labels = seq(0, max(d$y), by = .25), las = 1)
leg.text <- "Density over all tested animals"
legend(y = max(d$y) + .2, x = -.05, leg.text, lty = 1, bty = "n")
#dev.off()

ads <- subset(full.data, Age.Class == "Adult")
lambs <- subset(full.data, Age.Class == "Lamb")
yrs <- subset(full.data, Age.Class == "Yearling")

length(levels(factor(lambs$id)))
length(levels(factor(yrs$id)))

levels(factor(ads$id))[which(levels(factor(ads$id)) %in% levels(factor(lambs$id)))]
levels(factor(yrs$id))[which(levels(factor(yrs$id)) %in% levels(factor(lambs$id)))]
levels(factor(yrs$id))[which(levels(factor(yrs$id)) %in% levels(factor(ads$id)))]

hist(table(full.data$id), col = "grey80", main = "", xlab = "Number samples")
length(which(table(full.data$id)>=3))
length(which(table(full.data$id)>=5))

