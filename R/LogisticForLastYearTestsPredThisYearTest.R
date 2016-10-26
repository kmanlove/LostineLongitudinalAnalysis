# read in Lostine data
data <- read.csv("../Data/LostineMoviPastLung_150630.csv", header = T)

#-- build and explore determination cohort --#
Animal.tab <- table(data$Animal_ID)
AnimalByOutcome.tab <- table(data$Animal_ID, data$WADDL.Movi.PCR)
Animal.tab[which(Animal.tab >= 4)]
Animal.tab[which(Animal.tab == 3)]

ManySamplesOutcomes.tab <- AnimalByOutcome.tab[which(Animal.tab >= 4), ]
ManySamplesProps <- ManySamplesOutcomes.tab[ ,1] / (ManySamplesOutcomes.tab[ ,1] + ManySamplesOutcomes.tab[ ,3])

# histogram of prop positive
data$new.date <- strptime(data$Capture.Date, forma = "%m/%d/%Y")
plot(x = 0, y = 0, cex = 0, ylim = c(0, length(levels(factor(data$Animal_ID)))))
for(i in 1:length(levels(factor(data$Animal_ID)))){
  k <- subset(data, Animal_ID == levels(factor(data$Animal_ID))[i])
  points(y = rep(i, dim(k)[1]), x = k$new.date, pch = 16, col = as.numeric(k$WADDL.MoviPCR))
}

cols.in <- c("red", "pink", "grey50")
plot(data$new.date, as.numeric(data$Animal_ID), type = "p", xaxt = "n", pch = 16, col = cols.in[as.numeric(as.factor(data$WADDL.Movi.PCR))], xlab = "date", ylab = "Animal")
r <- as.POSIXct(round(range(data$new.date), "days"))
axis.Date(1, at = seq(r[1], r[2], by = "day"), format = "%H")

par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2), oma = c(1, 1, 1, 1))
hist(ManySamplesProps, breaks = 10, col = "grey80", xlab = "Proportion of PCR-positive tests for \n the 22 animals tested 4 or more times", main = "")
abline(v = 0.67, col = "red", lty = 2, lwd = 2)
text(x = 0.3, y = 4.5, "Consistently \n negative")
text(x = 0.8, y = 4.5, "Consistently \n positive")


# histograms of ELISA and qPCRs
par(mfrow = c(1, 2))
hist(data$Movipneumonia.ELISA, col = "grey80", xlab = "Titer", main = "Immune Response")
hist(data$MoviDirect_qPCR, col = "grey80", main = "Pathogen Load", xlab = "qPCR cycle threshold", breaks = 10)

# regression 1: logistic regression of first test outcome each year (pos / not) 
# by # positive tests last year | # tests last year

# 1. build data: rows are sheep-years
#   -- get number of sheep years
sheep.yrs <- unique(cbind(data$Animal_ID, data$CAP_BIOYR)) # = 134 in Lostine dataset

#   -- which sheep sampled in multiple years?
# mult.yr <- sheep.yrs[sheep.yrs[, 1] %in% which(table(sheep.yrs[, 1]) != 1), ]
# mult.yr[, 1] <- as.character(levels(factor(data$Animal_ID))[mult.yr[, 1]])
# 
# anim.id <- yr <- this.yr.result.1 <- last.yr.num.pos <- last.yr.num.tests <- vector("list", length(levels(factor(mult.yr[, 1]))))
# 
# for(i in 1:dim(mult.yr)[1]){
#   k <- subset(data, Animal_ID == mult.yr[i, 1])
#   k.yrs <- levels(factor(k$CAP_BIOYR))
#   # extract streaks of consecutive years
#   diff(as.numeric(as.character(k.yrs)))
#   yrs.to.use <- k.yrs[which(diff(as.numeric(as.character(k.yrs))) == 1)]
#   yr[[i]] <- this.yr.result.1[[i]] <- last.yr.num.pos[[i]] <- last.yr.num.tests[[i]] <- rep(NA, length(yrs.to.use))
#   # loop over years with following year's data available
#   for(j in 1:length(yrs.to.use)){
#     yr[[i]][j] <- yrs.to.use[j]
#     this.yr.result.1[[i]][j] <- ifelse(subset(k, CAP_BIOYR == (as.numeric(as.character(yrs.to.use))[j] + 1))$WADDL.Movi.PCR[1] == "Detected", 1, 0)
#     last.yr.num.tests[[i]][j] <- dim(subset(k, CAP_BIOYR == (yrs.to.use[j])))[1]
#     last.yr.num.pos[[i]][j] <- dim(subset(subset(k, CAP_BIOYR == (yrs.to.use[j])), WADDL.Movi.PCR == "Detected"))[1]
#   }
#   anim.id[[i]] <- rep(mult.yr[i, 1], length(yr[[i]]))
# }

anim.id <- yr <- this.yr.result.1 <- last.yr.num.pos <- last.yr.num.tests <- vector("list", length(levels(factor(data$Animal_ID))))

for(i in 1:length(levels(data$Animal_ID))){
  k <- subset(data, Animal_ID == levels(data$Animal_ID)[i])
  k.yrs <- levels(factor(k$CAP_BIOYR))
  # extract streaks of consecutive years
  diff(as.numeric(as.character(k.yrs)))
  yrs.to.use <- k.yrs[which(diff(as.numeric(as.character(k.yrs))) == 1)]
  yr[[i]] <- this.yr.result.1[[i]] <- last.yr.num.pos[[i]] <- last.yr.num.tests[[i]] <- rep(NA, length(yrs.to.use))
  # loop over years with following year's data available
  for(j in 1:length(yrs.to.use)){
    yr[[i]][j] <- yrs.to.use[j]
    this.yr.result.1[[i]][j] <- ifelse(subset(k, CAP_BIOYR == (as.numeric(as.character(yrs.to.use))[j] + 1))$WADDL.Movi.PCR[1] == "Detected", 1, 0)
    last.yr.num.tests[[i]][j] <- dim(subset(k, CAP_BIOYR == (yrs.to.use[j])))[1]
    last.yr.num.pos[[i]][j] <- dim(subset(subset(k, CAP_BIOYR == (yrs.to.use[j])), WADDL.Movi.PCR == "Detected"))[1]
  }
  anim.id[[i]] <- rep(levels(data$Animal_ID)[i], length(yr[[i]]))
}

anim.ids <- unlist(anim.id)
yrs <- unlist(yr)
this.yr.results <- unlist(this.yr.result.1)
last.yr.tests <- unlist(last.yr.num.tests)
last.yr.pos <- unlist(last.yr.num.pos)

last.yr.prop.pos <- last.yr.pos / last.yr.tests
yr.to.yr.dat <- as.data.frame(cbind(anim.ids, yrs, this.yr.results, last.yr.tests, last.yr.pos, last.yr.prop.pos))
yr.to.yr.dat <- subset(yr.to.yr.dat, is.na(yrs) == F)
yr.to.yr.dat$last.yr.prop.pos <- as.numeric(as.character(yr.to.yr.dat$last.yr.prop.pos))
yr.to.yr.dat$last.yr.pos <- as.numeric(as.character(yr.to.yr.dat$last.yr.pos))

# animals included:
levels(factor(yr.to.yr.dat$anim.id))

# 9 potential inclusion discrepancies 
#     (isolated year always excluded in this analysis):
# 06LO61 sampled in 2010, then not again until 2012
# 08LO29 sampled in 2009, then not again until 2012
# 10LO33 sampled in 2009, then not again until 2012
# 11LO40 sampled in 2010, then not again until 2012
# 11LO44 sampled in 2010, then not again until 2012
# 12LO27 sampled in 2012, then not again until 2014
# 13LO82 sampled in 2012, then not again until 2014
# 99LO03 sampled in 2009, then not again until 2012
# 99LO09 sampled in 2009, then not again until 2012
#     (three of these animals -- 12LO27, 13LO82, 99LO03 --
#     completely excluded from KRM analysis)
                                       
pos.last.yr.fit <- glm(this.yr.results ~ as.numeric(as.character(last.yr.pos)), 
                       offset = as.numeric(as.character(last.yr.tests)), 
                       family = "binomial", data = yr.to.yr.dat)

summary(pos.last.yr.fit)
yr.to.yr.dat$last.yr.prop.pos <- as.numeric(as.character(last.yr.prop.pos))

pos.last.yr.prop.fit <- glm(this.yr.results ~ last.yr.prop.pos, 
                       family = "binomial", data = yr.to.yr.dat)
summary(pos.last.yr.prop.fit)
pred.data <- as.data.frame(c(0, .5, 1))
names(pred.data) <- "last.yr.prop.pos"
prop.preds <- predict(pos.last.yr.prop.fit, newdata = pred.data, type = "response",
                      se.fit = T)
# prob of being pos on first test this year is 0.94 | pos on 100% of tests last year
odds.prop <- exp(2.9) 
odds.prop.int <- c(exp(2.9 - 1.96 * .795), exp(2.9 + 1.96 * .795))


pos.last.yr.count.fit <- glm(this.yr.results ~ last.yr.pos, 
                       family = "binomial", data = yr.to.yr.dat)
summary(pos.last.yr.count.fit)
pred.data.count <- as.data.frame(c(0, 1, 2, 3))
names(pred.data.count) <- "last.yr.pos"
prop.preds.count <- predict(pos.last.yr.count.fit, newdata = pred.data.count, 
                      type = "response", se.fit = T)


probs.count <- exp(1.10) / (1 + exp(1.10))
odds.count <- exp(1.10)
odds.count.int <- c(exp(1.10 - 1.96 * 0.4211), exp(1.10 + 1.96 * 0.4211))

# the odds an animal who is pos on 100% of tests last year being pos on first 
# test this year are 18.2 X the odds of an animal who is 100% neg last year being 
# pos on first test this year. 


pos.this.yr <- subset(yr.to.yr.dat, this.yr.results == "1")
pos.this.yr[which(pos.this.yr$last.yr.prop.pos < 0.2), ]
neg.this.yr <- subset(yr.to.yr.dat, this.yr.results == "0")
neg.this.yr[which(neg.this.yr$last.yr.prop.pos > 0.2), ]

par(mfrow = c(2, 1))
hist(pos.this.yr$last.yr.prop.pos, xlim = c(0, 1), breaks = 10, col = "grey80", 
     main = "Pos this year", xlab = "Prop pos last year")
hist(neg.this.yr$last.yr.prop.pos, xlim = c(0, 1), breaks = 10, col = "grey80",
     main = "Neg this year", xlab = "Prop post last year")