# code to build the prevalence by age color-block figure. 

full.data <- read.csv("./Data/LostineSample2016Final_WithClasses_20160927.csv", header = T)
full.data <- subset(full.data, movi_qpcr %in% c("Detected", "Not detected"))
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")
full.data <- subset(full.data, cap_bioyr >= 2008)

full.data$DaysSinceFirstCap <- full.data$NumberTests <- full.data$DaysSinceFirstOverallCap <- full.data$AgeRank <- full.data$leq5.ind <- rep(NA, dim(full.data)[1])
for(i in 1:dim(full.data)[1]){
  k <- subset(full.data, id == full.data$id[i])
  full.data$DaysSinceFirstCap[i] <- difftime(as.POSIXlt(strptime(full.data$capture_date[i], format = "%m/%d/%Y")),
                                             min(as.POSIXlt(strptime(k$capture_date, format = "%m/%d/%Y"))))
  full.data$DaysSinceFirstOverallCap[i] <- difftime(as.POSIXlt(strptime(full.data$capture_date[i], format = "%m/%d/%Y")),
                                             min(as.POSIXlt(strptime(full.data$capture_date, format = "%m/%d/%Y"))))
  full.data$NumberTests[i] <- dim(k)[1]
  full.data$AgeRank[i] <- min(k$age_capture)
  full.data$leq5.ind[i] <- ifelse(full.data$age_capture[i] <= 5, 1, 0)
}

full.data$pt.col <- ifelse(full.data$qPCRResult == "P", rgb(1, 0, 0, alpha = .7), 
                 ifelse(full.data$qPCRResult == "N", rgb(0, 0, 1, alpha = .7), 
                        rgb(.6, .6, .6, alpha = .5)))

reduced.data <- subset(full.data, NumberTests >= 2)

par(mfrow = c(1, 1), las = 1)
plot(as.numeric(factor(reduced.data$id)) ~ 
     as.numeric(as.character(reduced.data$DaysSinceFirstCap)),
     xlab = "Days since first capture", ylab = "Animal ID",
     main = "ELISAs through time within animals",
     cex = (reduced.data$movi_elisa/100 * 2), col = reduced.data$pt.col, pch = 16)
for(i in 1:length(levels(factor(reduced.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

plot(as.numeric(factor(reduced.data$id)) ~ 
       as.numeric(as.character(reduced.data$DaysSinceFirstOverallCap)),
     xlab = "Days since first capture", ylab = "Animal ID",
     cex = (reduced.data$movi_elisa/100 * 2), col = reduced.data$pt.col, pch = 16)
for(i in 1:length(levels(factor(reduced.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

plot(as.numeric(factor(reduced.data$id)) ~ 
       jitter(as.numeric(factor(reduced.data$qPCRResult)), .5),
     cex = (reduced.data$Movi.ELISA/100 * 2), col = reduced.data$pt.col, 
     pch = 16,
     xlim = c(0.5, 2.5), xaxt = "n",
     xlab = "", ylab = "Animal ID")
axis(side = 1, at = c(1, 2), labels = c("Negative", "Positive"))
for(i in 1:length(levels(factor(reduced.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

boxplot(reduced.data$movi_elisa ~ reduced.data$qPCRResult, col = "grey80",
        xaxt = "n", ylab = "ELISA titer")
axis(side = 1, at = c(1, 2), labels = c("qPCR-Negative", "qPCR-Positive"))
kruskal.ELISAByqPCR <- kruskal.test(as.numeric(as.character(reduced.data$movi_elisa)) ~ factor(reduced.data$qPCRResult))
t.test.ELISAbyqPCR <- t.test(as.numeric(as.character(reduced.data$movi_elisa)) ~ factor(reduced.data$qPCRResult))

plot(as.numeric(as.character(reduced.data$age_capture)) ~ 
       jitter(as.numeric(factor(reduced.data$qPCRResult)), .5),
     cex = (reduced.data$movi_elisa/100 * 2), col = reduced.data$pt.col, 
     pch = 16,
     xlim = c(0.5, 2.5), xaxt = "n",
     xlab = "", ylab = "Age at capture")
axis(side = 1, at = c(1, 2), labels = c("qPCR-Negative", "qPCR-Positive"))
abline(h = 5, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
abline(h = 10, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
abline(h = 15, lty = 2, lwd = .5, col = rgb(.6, .6, .6))

par(mfrow = c(1, 2))
plot(as.numeric(as.character(reduced.data$age_capture)) ~ 
       jitter(as.numeric(factor(reduced.data$qPCRResult)), .5),
     cex = (1/reduced.data$Direct.swab.qPCR * 40), col = reduced.data$pt.col, 
     pch = 16,
     main = "Direct swab qPCR",
     xlim = c(0.5, 2.5), xaxt = "n",
     xlab = "", ylab = "Age at capture")
axis(side = 1, at = c(1, 2), labels = c("qPCR-Negative", "qPCR-Positive"))
abline(h = 5, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
abline(h = 10, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
abline(h = 15, lty = 2, lwd = .5, col = rgb(.6, .6, .6))

plot(as.numeric(as.character(reduced.data$age_capture)) ~ 
       jitter(as.numeric(factor(reduced.data$qPCRResult)), .5),
     cex = (1/reduced.data$Nasal.wash.qPCR * 40), col = reduced.data$pt.col, 
     pch = 16,
     main = "Nasal wash qPCR",
     xlim = c(0.5, 2.5), xaxt = "n",
     xlab = "", ylab = "Age at capture")
axis(side = 1, at = c(1, 2), labels = c("qPCR-Negative", "qPCR-Positive"))
abline(h = 5, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
abline(h = 10, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
abline(h = 15, lty = 2, lwd = .5, col = rgb(.6, .6, .6))

# par(mfrow = c(1, 1))
# plot(as.numeric(as.character(full.data$Direct.swab.qPCR)) ~ full.data$Nasal.wash.qPCR,
#      pch = ifelse(full.data$leq5.ind == 1, 16, 1),
#      xlab = "Nasal wash qPCR", ylab = "Direct swab qPCR", cex = 2)
# leg.text <- c("5 or younger", "Older than 5")
# legend("topleft", leg.text, pch = c(16, 1), bty = "n")


# age-spec prevalence
prev.tab <- prop.table(table(reduced.data$MostLikelyAge, reduced.data$pt.col), margin = 1)
prev.col <- rgb(prev.tab[,2], 0, prev.tab[,1], alpha = 1)

# svg("./Plots/AgeByPrev_IndividAndMarginal_20160507_V2.svg", height = 12, width = 10)
layout(matrix(c(1, 1, 1, 1, 2), ncol = 1, byrow = T))
plot(as.numeric(factor(reduced.data$id)) ~ 
       as.numeric(as.character(reduced.data$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal #",
     cex = 2, #cex = (reduced.data$Movi.ELISA/100 * 5), 
     col = reduced.data$pt.col, 
     pch = 16, xlim = c(0, 18), ylim = c(1, 45))
for(i in 1:length(levels(factor(reduced.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}
axis(side = 2, at = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45), 
     labels = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45), cex = .7)

plot(rep(1, 18) ~ seq(1,18, by = 1), col = prev.col, pch = 15, cex = 8, yaxt = "n",
     xlab = "Most likely age", ylab = "", xlim = c(0, 18))
# dev.off()

# age-spec prevalence
full.data$MostLikelyAge <- ifelse(full.data$MostLikelyAge == 0, 1, full.data$MostLikelyAge)
full.data$id <- ifelse(full.data$id == "99L03", "000L03",
                       ifelse(full.data$id == "99L09", "000L09",
                              as.character(full.data$id)))
prev.tab <- prop.table(table(full.data$MostLikelyAge, full.data$pt.col), margin = 1)
prev.col <- rgb(prev.tab[,2], 0, prev.tab[,1], alpha = 1)
prev.col.in <- c(prev.col[1:5], "white", prev.col[6:18])

reduced.data$AgeOrder <- rank(reduced.data$MostLikelyAge)
# sort animals by min(AgeatCapture)
min.ageatcapture <- rep(NA, length(levels(factor(full.data$id))))
for(i in 1:length(levels(factor(full.data$id)))){
  k <- subset(full.data, id == levels(factor(full.data$id))[i])
  min.ageatcapture[i] <- min(k$AgeatCapture)
}

order(min.ageatcapture)
cbind(min.ageatcapture, rank(min.ageatcapture, ties = "first"))
full.data$animal.order <- rank(min.ageatcapture, ties = "first")[as.numeric(factor(full.data$id))]



layout(matrix(c(1, 1, 1, 1, 2), ncol = 1, byrow = T))
plot(rev(full.data$animal.order) ~ 
       as.numeric(as.character(full.data$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal #",
     cex = 2, #cex = (reduced.data$Movi.ELISA/100 * 5), 
     col = full.data$pt.col, 
     pch = 16, xlim = c(0, 18), ylim = c(1, 80))
for(i in 1:length(levels(factor(full.data$animal.order)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

plot(rep(1, 18) ~ seq(1, 18, by = 1), col = prev.col.in[1:18], pch = 15, cex = 8, yaxt = "n",
     xlab = "Most likely age", ylab = "", xlim = c(0, 18))
# dev.off()

#--------------------------------#
#-- COLOR BLOCK AGE PREVALENCE --#
#--------------------------------#

#svg("./Plots/AgeByPrev_IndividAndMarginal_20160929.svg", height = 7, width = 5)
layout(matrix(c(1, 1, 1, 1, 2), ncol = 1, byrow = T))
plot((as.numeric(factor(full.data$id))) ~ 
       as.numeric(as.character(full.data$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal #",
     cex = 2, #cex = (reduced.data$Movi.ELISA/100 * 5), 
     col = full.data$pt.col, 
     pch = 16, xlim = c(0, 18), ylim = c(1, 83))
for(i in 1:length(levels(factor(full.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

plot(rep(1, 18) ~ seq(1, 18, by = 1), col = prev.col.in[1:18], pch = 15, cex = 8, yaxt = "n",
     xlab = "Most likely age", ylab = "", xlim = c(0, 18))
#dev.off()
#------------------------------------#
#-- END COLOR BLOCK AGE-PREVALENCE --#
#------------------------------------#

# ELISA version
# age-spec prevalence
reduced.data$ELISA.ind <- ifelse(reduced.data$Movi.ELISA >= 50, 2, 1)
reduced.data.2 <- subset(reduced.data, is.na(ELISA.tab) == F)
ELISA.tab <- tapply(reduced.data.2$Movi.ELISA, reduced.data.2$MostLikelyAge, mean)
ELISA.col <- rep(NA, length(ELISA.tab))
for(i in 1:length(ELISA.tab)){
  ELISA.col[i] <- ifelse(is.na(ELISA.tab[i]) == T, "grey60", rgb(ELISA.tab[i]/100, 0, (1 - ELISA.tab[i]/100), alpha = 1))
}

# svg("./Plots/AgeByELISA_IndividAndMarginal_20160505.svg", height = 12, width = 10)
layout(matrix(c(1, 1, 1, 1, 2), ncol = 1, byrow = T))
plot(as.numeric(factor(reduced.data$id)) ~ 
       as.numeric(as.character(reduced.data$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal ID",
     cex = (reduced.data$Movi.ELISA/100 * 2), col = reduced.data$ELISA.ind, 
     pch = 16, xlim = c(0, 18))
for(i in 1:length(levels(factor(reduced.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

reduced.data$ELISA.col.ind <- rep(NA, dim(reduced.data)[1])
for(i in 1:dim(reduced.data)[1]){
  reduced.data$ELISA.col.ind[i] <- ifelse(is.na(reduced.data$Movi.ELISA[i]) == T, "grey60", 
                                          rgb(reduced.data$Movi.ELISA[i] / 100, 0, (1 - reduced.data$Movi.ELISA[i] / 100), alpha = 1))
}

# ELISA Full data
full.data$ELISA.ind <- ifelse(full.data$Movi.ELISA >= 50, 2, 1)
full.data.2 <- subset(full.data, is.na(ELISA.tab) == F)
ELISA.tab <- tapply(full.data.2$Movi.ELISA, full.data.2$MostLikelyAge, mean)
ELISA.col <- rep(NA, length(ELISA.tab))
for(i in 1:length(ELISA.tab)){
  ELISA.col[i] <- ifelse(is.na(ELISA.tab[i]) == T, "grey60", 
                         rgb(ELISA.tab[i]/100, 0, (1 - ELISA.tab[i]/100), alpha = 1))
}

# svg("./Plots/AgeByELISA_IndividAndMarginal_20160505.svg", height = 12, width = 10)
layout(matrix(c(1, 1, 1, 1, 2), ncol = 1, byrow = T))
plot(as.numeric(factor(full.data$id)) ~ 
       as.numeric(as.character(full.data$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal ID",
     cex = (full.data$Movi.ELISA/100 * 5), col = full.data$ELISA.ind, 
     pch = 16, xlim = c(0, 18))
for(i in 1:length(levels(factor(full.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

full.data$ELISA.col.ind <- rep(NA, dim(full.data)[1])
for(i in 1:dim(full.data)[1]){
  full.data$ELISA.col.ind[i] <- ifelse(is.na(full.data$Movi.ELISA[i]) == T, "grey60", 
#                                       rgb(full.data$Movi.ELISA[i] / 100, 0, (1 - full.data$Movi.ELISA[i] / 100), alpha = 1)
                                       rgb(full.data$Movi.ELISA[i] / 100, 0, 0, alpha = 1)
  )
}


# ELISA version
# age-spec prevalence
full.data$ELISA.ind <- ifelse(reduced.data$Movi.ELISA >= 50, 2, 1)
reduced.data.2 <- subset(reduced.data, is.na(ELISA.tab) == F)
ELISA.tab <- tapply(reduced.data.2$Movi.ELISA, reduced.data.2$MostLikelyAge, mean)

full.elisa <- subset(full.data, is.na(Movi.ELISA) == F)

plot(Movi.ELISA ~ MostLikelyAge, data = full.elisa, xlim = c(0, 20),
     xlab = "Animal age", ylab = "Movi ELISA (% Inhibition)", pch = 16, col = rgb(0, 0, 0, alpha = .4))
lines(lowess(full.elisa$Movi.ELISA ~ full.elisa$MostLikelyAge), lwd = 1.5)

ELISA.col <- rep(NA, length(ELISA.tab))
for(i in 1:length(ELISA.tab)){
  ELISA.col[i] <- ifelse(is.na(ELISA.tab[i]) == T, "grey60", rgb(ELISA.tab[i]/100, 0, (1 - ELISA.tab[i]/100), alpha = 1))
}

plot(rep(1, 18) ~ seq(1,18, by = 1), col = ELISA.col, pch = 15, cex = 8, yaxt = "n",
     xlab = "Most likely age", ylab = "", xlim = c(0, 18))

# svg("./Plots/AgeByELISA_IndividAndMarginal_20160505.svg", height = 12, width = 10)
layout(matrix(c(1, 1, 1, 1, 2), ncol = 1, byrow = T))
plot(as.numeric(factor(reduced.data$id)) ~ 
       as.numeric(as.character(reduced.data$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal ID",
     cex = (reduced.data$Movi.ELISA/100 * 5), col = reduced.data$ELISA.ind, 
     pch = 16, xlim = c(0, 18))
for(i in 1:length(levels(factor(reduced.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

reduced.data$ELISA.col.ind <- rep(NA, dim(reduced.data)[1])
for(i in 1:dim(reduced.data)[1]){
  reduced.data$ELISA.col.ind[i] <- ifelse(is.na(reduced.data$Movi.ELISA[i]) == T, "grey60", 
                                          rgb(reduced.data$Movi.ELISA[i] / 100, 0, (1 - reduced.data$Movi.ELISA[i] / 100), alpha = 1))
}

plot(rep(1, 18) ~ seq(1,18, by = 1), col = ELISA.col, pch = 15, cex = 8, yaxt = "n",
     xlab = "Most likely age", ylab = "", xlim = c(0, 18))
# dev.off()


layout(matrix(c(1, 1, 1, 1, 2), ncol = 1, byrow = T))
plot(as.numeric(factor(reduced.data$id)) ~ 
       as.numeric(as.character(reduced.data$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal ID",
     cex = (reduced.data$Movi.ELISA/100 * 5), col = reduced.data$ELISA.col.ind, 
     pch = 16, xlim = c(0, 18), ylim = c(1, 44))
for(i in 1:length(levels(factor(reduced.data$id)))){
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
}

plot(rep(1, 18) ~ seq(1,18, by = 1), col = ELISA.col, pch = 15, cex = 8, yaxt = "n",
     xlab = "Most likely age", ylab = "", xlim = c(0, 18))


lambs.yrs.only <- subset(full.data, AgeatEntry <= 1.9 & NumberTests >= 2 & ENTRY_BIOYR >= 2010)
plot(as.numeric(factor(lambs.yrs.only$Animal.ID)) ~ 
       as.numeric(as.character(lambs.yrs.only$MostLikelyAge)),
     xlab = "Most likely age", ylab = "Animal ID",
     cex = (lambs.yrs.only$Movi.ELISA/100 * 5), col = lambs.yrs.only$pt.col, 
     pch = 16, xlim = c(0, 3))
for(i in 1:length(levels(factor(lambs.yrs.only$Animal.ID)))){
  k <- subset(full.data, as.character(Animal.ID) == levels(factor(lambs.yrs.only$Animal.ID))[i])
  abline(h = i, lty = 2, lwd = .5, col = rgb(.6, .6, .6))
  segments(x0 = min(k$MostLikelyAge), x1 = max(k$MostLikelyAge), y0 = i, y1 = i)
}
