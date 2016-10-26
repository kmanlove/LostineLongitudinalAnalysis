movi <- read.csv("./Data/MoviSurvivalDate_20160220.csv", header = T)

movi$stopdate.Date <- as.Date(movi$stopdate, format = "%m/%d/%Y")
movi$startdate.Date <- as.Date(movi$startdate, format = "%m/%d/%Y")
movi$previousstartsamp.Date <- as.Date(movi$previousstartsamp, format = "%m/%d/%Y")
movi$previousendsamp.Date <- as.Date(movi$previousendsamp, format = "%m/%d/%Y")

movi$start.mid <- movi$previousstartsamp.Date + floor(movi$startdate.Date - movi$previousstartsamp.Date)/2
movi$start.mid[is.na(movi$previousstartsamp.Date) == T] <- movi$startdate.Date[is.na(movi$previousstartsamp.Date) == T]
movi$stop.mid <- movi$previousendsamp.Date + floor(movi$stopdate.Date - movi$previousendsamp.Date)/2
movi$stop.mid[is.na(movi$previousendsamp.Date) == T] <- movi$stopdate.Date[is.na(movi$previousendsamp.Date) == T]
movi$diffmid <- difftime(movi$stop.mid, movi$start.mid)


timepos <- subset(movi, origstatus == "pos")
timepos$event1 <- ifelse(timepos$endstatus == "pos", 0, 1)
require(survival)

testfit <- survfit(Surv(rep(0, dim(timepos)[1]), as.numeric(as.character(timepos$diffmid)), 
                        event = timepos$event1) ~ 1)
plot(testfit)

noind <- subset(movi, origstatus != "ind" & endstatus != "ind")

testfit2 <- survfit(Surv(rep(0, dim(noind)[1]), as.numeric(as.character(noind$diffmid)), 
                        event = noind$rightcensor) ~ factor(noind$origstatus))
plot(testfit2, conf.int = T, col = c("red", "blue"))

timepos <- subset(noind, origstatus == "pos")
posfit <- survfit(Surv(rep(0, dim(timepos)[1]), as.numeric(as.character(timepos$diffmid)), 
                         event = timepos$rightcensor) ~ 1)

timeneg <- subset(noind, origstatus == "neg")
negfit <- survfit(Surv(rep(0, dim(timeneg)[1]), as.numeric(as.character(timeneg$diffmid)), 
                       event = timeneg$rightcensor) ~ 1)

par(mfrow = c(1, 2), las = 1)
plot(posfit, conf.int = T, col = c("red"), xlim = c(0, 1000),
     xlab = "Elapsed time since first test (days)", ylab = "Proportion still positive")
abline(v = 365, lty = 3, col = "grey40")
abline(v = 730, lty = 3, col = "grey40")
text(x = 850, y = .95, "N = 27")
text(x = 365, y = .01, cex = .6, "1 year")
text(x = 730, y = .01, cex = .6, "2 years")

plot(negfit, conf.int = T, col = c("blue"), xlim = c(0, 1000),
     xlab = "Elapsed time since first test (days)", ylab = "Proportion still negative")
abline(v = 365, lty = 3, col = "grey40")
abline(v = 730, lty = 3, col = "grey40")
text(x = 850, y = .95, "N = 19")
text(x = 365, y = .01, cex = .6, "1 year")
text(x = 730, y = .01, cex = .6, "2 years")




