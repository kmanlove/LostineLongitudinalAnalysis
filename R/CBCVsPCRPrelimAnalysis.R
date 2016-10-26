full.data <- read.csv("./Data/LostineSampleData_20160505_v1.csv", header = T)
full.data$id <- full.data$Animal.ID
full.data$movi_qpcr <- full.data$WADDL.Movi.PCR
full.data$cap_bioyr <- full.data$CAP_BIOYR
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")
full.data <- subset(full.data, CAP_BIOYR >= 2008)

full.data$DaysSinceFirstCap <- full.data$NumberTests <- full.data$DaysSinceFirstOverallCap <- full.data$AgeRank <- full.data$leq5.ind <- rep(NA, dim(full.data)[1])
for(i in 1:dim(full.data)[1]){
  k <- subset(full.data, Animal.ID == full.data$Animal.ID[i])
  full.data$DaysSinceFirstCap[i] <- difftime(as.POSIXlt(strptime(full.data$Capture.Date[i], format = "%d/%m/%Y")),
                                             min(as.POSIXlt(strptime(k$Capture.Date, format = "%d/%m/%Y"))))
  full.data$DaysSinceFirstOverallCap[i] <- difftime(as.POSIXlt(strptime(full.data$Capture.Date[i], format = "%d/%m/%Y")),
                                                    min(as.POSIXlt(strptime(full.data$Capture.Date, format = "%d/%m/%Y"))))
  full.data$NumberTests[i] <- dim(k)[1]
  full.data$AgeRank[i] <- min(k$AgeatCapture)
  full.data$leq5.ind[i] <- ifelse(full.data$AgeatCapture[i] <= 5, 1, 0)
  print(i)
}



full.data$pt.col <- ifelse(full.data$qPCRResult == "P", rgb(1, 0, 0, alpha = .7), 
                           ifelse(full.data$qPCRResult == "N", rgb(0, 0, 1, alpha = .7), 
                                  rgb(.6, .6, .6, alpha = .5)))

full.data$Lung.Larvae..lpg. <- ifelse(full.data$Lung.Larvae..lpg. == "Negative" | full.data$Lung.Larvae..lpg. == "negative", 0,
                               ifelse(full.data$Lung.Larvae..lpg. == "insufficient sample" | full.data$Lung.Larvae..lpg. == "insuficient sample", "",
                                      full.data$Lung.Larvae..lpg.))
worms.small <- subset(full.data, !(Lung.Larvae..lpg.== "") & is.na(full.data$qPCRResult) == FALSE)
worms.small.juv <- subset(worms.small, MostLikelyAge <= 2)
worms.small.ad <- subset(worms.small, MostLikelyAge > 2)

par(mfrow = c(1, 2), cex.axis = .8, cex.lab = 1)
boxplot((as.numeric(as.character(worms.small.juv$Lung.Larvae..lpg.)) + 1) ~ as.numeric(factor(worms.small.juv$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), log = "y",
        main = "Animals <= 2", ylab = "Worm burden")
boxplot((as.numeric(as.character(worms.small.ad$Lung.Larvae..lpg.)) + 1) ~ as.numeric(factor(worms.small.ad$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), log = "y", main = "Animals older than 2",
        ylab = "Worm burden")

blood.small <- subset(full.data, is.na(LymphPct) == F & is.na(full.data$qPCRResult) == FALSE)
blood.small.juv <- subset(blood.small, MostLikelyAge <= 2)
blood.small.ad <- subset(blood.small, MostLikelyAge > 2)

par(mfrow = c(3, 4), cex.axis = .8, cex.lab = 1)
# Band Pct
boxplot((as.numeric(as.character(blood.small.juv$BandPct)) + 1) ~ as.numeric(factor(blood.small.juv$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), 
        main = "Animals <= 2", ylab = "BandPct")
boxplot((as.numeric(as.character(blood.small.ad$BandPct)) + 1) ~ as.numeric(factor(blood.small.ad$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), main = "Animals older than 2",
        ylab = "BandPct")
# Seg Pct
boxplot((as.numeric(as.character(blood.small.juv$SegPct)) + 1) ~ as.numeric(factor(blood.small.juv$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), 
        main = "Animals <= 2", ylab = "SegPct")
boxplot((as.numeric(as.character(blood.small.ad$SegPct)) + 1) ~ as.numeric(factor(blood.small.ad$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), main = "Animals older than 2",
        ylab = "SegPct")
# Lymph Pct
boxplot((as.numeric(as.character(blood.small.juv$LymphPct)) + 1) ~ as.numeric(factor(blood.small.juv$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), 
        main = "Animals <= 2", ylab = "LymphPct")
boxplot((as.numeric(as.character(blood.small.ad$LymphPct)) + 1) ~ as.numeric(factor(blood.small.ad$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), main = "Animals older than 2",
        ylab = "LymphPct")
# Mono Pct
boxplot((as.numeric(as.character(blood.small.juv$MonoPct)) + 1) ~ as.numeric(factor(blood.small.juv$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), 
        main = "Animals <= 2", ylab = "MonoPct")
boxplot((as.numeric(as.character(blood.small.ad$MonoPct)) + 1) ~ as.numeric(factor(blood.small.ad$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), main = "Animals older than 2",
        ylab = "MonoPct")
# Eosin Pct
boxplot((as.numeric(as.character(blood.small.juv$EosinPct)) + 1) ~ as.numeric(factor(blood.small.juv$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), 
        main = "Animals <= 2", ylab = "EosinPct")
boxplot((as.numeric(as.character(blood.small.ad$EosinPct)) + 1) ~ as.numeric(factor(blood.small.ad$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), main = "Animals older than 2",
        ylab = "EosinPct")
# Basophil Pct
boxplot((as.numeric(as.character(blood.small.juv$BasophilPct)) + 1) ~ as.numeric(factor(blood.small.juv$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), 
        main = "Animals <= 2", ylab = "BasophilPct")
boxplot((as.numeric(as.character(blood.small.ad$BasophilPct)) + 1) ~ as.numeric(factor(blood.small.ad$qPCRResult)),
        col = "grey60", names = c("Movi-neg", "Movi-pos"), main = "Animals older than 2",
        ylab = "BasophilPct")

#-----------------------------#
#--- ELISA vs. PCR results ---#
#-----------------------------#

plot(as.numeric(as.character(full.data$Movi.ELISA)) ~ jitter(as.numeric(factor(full.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = "black", xaxt = "n", xlab = "Movi PCR", ylab = "Movi ELISA (% Inhibition)", las = 1, xlim = c(0.5, 3.5))
abline(h = 40, lty = 2)
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))

full.data$ELISA.ind <- ifelse(full.data$Movi.ELISA <= 40, 0, 1)

chisq.test(x = full.data$ELISA.ind, y = full.data$WADDL.Movi.PCR,
           simulate.p.value = T, B = 2000)
chisq.test(x = full.data$ELISA.ind, y = full.data$qPCRResult,
           simulate.p.value = T, B = 2000)
chisq.test(x = full.data$ELISA.ind, y = full.data$qPCRResult,
           simulate.p.value = F)

