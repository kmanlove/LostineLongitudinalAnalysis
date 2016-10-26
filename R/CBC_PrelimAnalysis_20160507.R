cbc.data <- read.csv("./Data/BHS_CBC_Data.csv", sep = ",")


cbc.small <- subset(cbc.data, select = c("Herd", "Sex", "ACLASS_SAMP", "WADDL.Movi.PCR", 
                                         "Direct.swab.qPCR", "Movipneumonia.ELISA",
                                         "WBCx1000peruL", "SegPct", 
                                         "LymphPct", "MonoPct", "EosinPct", 
                                         "Fibrinogen", "Hgb.gperdL",
                                         "PCVPct", "MCVfl", "RDW_Pct", "Plateletx1000peruL"))

# now, cut out Direct.swab.qPCR.
cbc.small <- subset(cbc.data, select = c("Herd", "Sex", "ACLASS_SAMP", "WADDL.Movi.PCR", 
                                         "Movipneumonia.ELISA",
                                         "WBCx1000peruL", "SegPct", 
                                         "LymphPct", "MonoPct", "EosinPct", 
                                         "Fibrinogen", "Hgb.gperdL",
                                         "PCVPct", "MCVfl", "RDW_Pct", "Plateletx1000peruL"))


cbc.small.complete.cases <- complete.cases(cbc.small)
cbc.small.complete <- cbc.small[cbc.small.complete.cases, ]
cbc.small.complete$qPCRResult <- ifelse(cbc.small.complete$WADDL.Movi.PCR == "Detected" | cbc.small.complete$WADDL.Movi.PCR == "Indeterminate", 2, 1)

pca.out <- prcomp(cbc.small.complete[, -c(1:5, 17)], scale. = T, center = T)

plot(pca.out$x[, 1] ~ pca.out$x[, 2], col = cbc.small.complete$qPCRResult, pch = 16,
     cex = cbc.small.complete$Movipneumonia.ELISA / 15)



round(pca.out$rotation, 2)
cbc.small.complete$PC1 <- pca.out$x[, 1]
cbc.small.complete$PC2 <- pca.out$x[, 2]

# Plot PC loadings
par(mfrow = c(2, 2), las = 1, cex.axis = .8, oma = c(0, 2, 0, 0), mar = c(6, 6, 1, 1))
plot(1 ~ 1, cex = 0, xlim = c(-.6, .6), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC1", bty = "n")
polygon(x = c(rep(pca.out$rotation[1, 1], 2), rep(0, 2)), y = c(11, 10, 10, 11), col = "grey80")
polygon(x = c(rep(pca.out$rotation[2, 1], 2), rep(0, 2)), y = c(10, 9, 9, 10), col = "grey80")
polygon(x = c(rep(pca.out$rotation[3, 1], 2), rep(0, 2)), y = c(9, 8, 8, 9), col = "grey80")
polygon(x = c(rep(pca.out$rotation[4, 1], 2), rep(0, 2)), y = c(8, 7, 7, 8), col = "grey80")
polygon(x = c(rep(pca.out$rotation[5, 1], 2), rep(0, 2)), y = c(7, 6, 6, 7), col = "grey80")
polygon(x = c(rep(pca.out$rotation[6, 1], 2), rep(0, 2)), y = c(6, 5, 5, 6), col = "grey80")
polygon(x = c(rep(pca.out$rotation[7, 1], 2), rep(0, 2)), y = c(5, 4, 4, 5), col = "grey80")
polygon(x = c(rep(pca.out$rotation[8, 1], 2), rep(0, 2)), y = c(4, 3, 3, 4), col = "grey80")
polygon(x = c(rep(pca.out$rotation[9, 1], 2), rep(0, 2)), y = c(3, 2, 2, 3), col = "grey80")
polygon(x = c(rep(pca.out$rotation[10, 1], 2), rep(0, 2)), y = c(2, 1, 1, 2), col = "grey80")
polygon(x = c(rep(pca.out$rotation[11, 1], 2), rep(0, 2)), y = c(1, 0, 0, 1), col = "grey80")
abline(v = 0, lwd = 3)
axis(side = 2, at = c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), c("Plateletx1000perul", "RDW_Pct", "MCVfl", "PCVPct","Hgb.gperdL", "Fibrinogen", "EosinPct","MonoPct", "LymphPct","SegPct", "WBCx1000peruL"))
abline(v = .25, lty = 2, col = "grey20")
abline(v = -.25, lty = 2, col = "grey20")

plot(1 ~ 1, cex = 0, xlim = c(-.6, .6), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC2", bty = "n")
polygon(x = c(rep(pca.out$rotation[1, 2], 2), rep(0, 2)), y = c(11, 10, 10, 11), col = "grey80")
polygon(x = c(rep(pca.out$rotation[2, 2], 2), rep(0, 2)), y = c(10, 9, 9, 10), col = "grey80")
polygon(x = c(rep(pca.out$rotation[3, 2], 2), rep(0, 2)), y = c(9, 8, 8, 9), col = "grey80")
polygon(x = c(rep(pca.out$rotation[4, 2], 2), rep(0, 2)), y = c(8, 7, 7, 8), col = "grey80")
polygon(x = c(rep(pca.out$rotation[5, 2], 2), rep(0, 2)), y = c(7, 6, 6, 7), col = "grey80")
polygon(x = c(rep(pca.out$rotation[6, 2], 2), rep(0, 2)), y = c(6, 5, 5, 6), col = "grey80")
polygon(x = c(rep(pca.out$rotation[7, 2], 2), rep(0, 2)), y = c(5, 4, 4, 5), col = "grey80")
polygon(x = c(rep(pca.out$rotation[8, 2], 2), rep(0, 2)), y = c(4, 3, 3, 4), col = "grey80")
polygon(x = c(rep(pca.out$rotation[9, 2], 2), rep(0, 2)), y = c(3, 2, 2, 3), col = "grey80")
polygon(x = c(rep(pca.out$rotation[10, 2], 2), rep(0, 2)), y = c(2, 1, 1, 2), col = "grey80")
polygon(x = c(rep(pca.out$rotation[11, 2], 2), rep(0, 2)), y = c(1, 0, 0, 1), col = "grey80")
abline(v = 0, lwd = 3)
axis(side = 2, at = c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), c("Plateletx1000perul", "RDW_Pct", "MCVfl", "PCVPct","Hgb.gperdL", "Fibrinogen", "EosinPct","MonoPct", "LymphPct","SegPct", "WBCx1000peruL"))
abline(v = .25, lty = 2, col = "grey20")
abline(v = -.25, lty = 2, col = "grey20")


plot(1 ~ 1, cex = 0, xlim = c(-.6, .6), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC1", bty = "n")
polygon(x = c(rep(pca.out$rotation[1, 3], 2), rep(0, 2)), y = c(11, 10, 10, 11), col = "grey80")
polygon(x = c(rep(pca.out$rotation[2, 3], 2), rep(0, 2)), y = c(10, 9, 9, 10), col = "grey80")
polygon(x = c(rep(pca.out$rotation[3, 3], 2), rep(0, 2)), y = c(9, 8, 8, 9), col = "grey80")
polygon(x = c(rep(pca.out$rotation[4, 3], 2), rep(0, 2)), y = c(8, 7, 7, 8), col = "grey80")
polygon(x = c(rep(pca.out$rotation[5, 3], 2), rep(0, 2)), y = c(7, 6, 6, 7), col = "grey80")
polygon(x = c(rep(pca.out$rotation[6, 3], 2), rep(0, 2)), y = c(6, 5, 5, 6), col = "grey80")
polygon(x = c(rep(pca.out$rotation[7, 3], 2), rep(0, 2)), y = c(5, 4, 4, 5), col = "grey80")
polygon(x = c(rep(pca.out$rotation[8, 3], 2), rep(0, 2)), y = c(4, 3, 3, 4), col = "grey80")
polygon(x = c(rep(pca.out$rotation[9, 3], 2), rep(0, 2)), y = c(3, 2, 2, 3), col = "grey80")
polygon(x = c(rep(pca.out$rotation[10, 3], 2), rep(0, 2)), y = c(2, 1, 1, 2), col = "grey80")
polygon(x = c(rep(pca.out$rotation[11, 3], 2), rep(0, 2)), y = c(1, 0, 0, 1), col = "grey80")
abline(v = 0, lwd = 3)
axis(side = 2, at = c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), c("Plateletx1000perul", "RDW_Pct", "MCVfl", "PCVPct","Hgb.gperdL", "Fibrinogen", "EosinPct","MonoPct", "LymphPct","SegPct", "WBCx1000peruL"))
abline(v = .25, lty = 2, col = "grey20")
abline(v = -.25, lty = 2, col = "grey20")

plot(1 ~ 1, cex = 0, xlim = c(-.6, .6), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC2", bty = "n")
polygon(x = c(rep(pca.out$rotation[1, 4], 2), rep(0, 2)), y = c(11, 10, 10, 11), col = "grey80")
polygon(x = c(rep(pca.out$rotation[2, 4], 2), rep(0, 2)), y = c(10, 9, 9, 10), col = "grey80")
polygon(x = c(rep(pca.out$rotation[3, 4], 2), rep(0, 2)), y = c(9, 8, 8, 9), col = "grey80")
polygon(x = c(rep(pca.out$rotation[4, 4], 2), rep(0, 2)), y = c(8, 7, 7, 8), col = "grey80")
polygon(x = c(rep(pca.out$rotation[5, 4], 2), rep(0, 2)), y = c(7, 6, 6, 7), col = "grey80")
polygon(x = c(rep(pca.out$rotation[6, 4], 2), rep(0, 2)), y = c(6, 5, 5, 6), col = "grey80")
polygon(x = c(rep(pca.out$rotation[7, 4], 2), rep(0, 2)), y = c(5, 4, 4, 5), col = "grey80")
polygon(x = c(rep(pca.out$rotation[8, 4], 2), rep(0, 2)), y = c(4, 3, 3, 4), col = "grey80")
polygon(x = c(rep(pca.out$rotation[9, 4], 2), rep(0, 2)), y = c(3, 2, 2, 3), col = "grey80")
polygon(x = c(rep(pca.out$rotation[10, 4], 2), rep(0, 2)), y = c(2, 1, 1, 2), col = "grey80")
polygon(x = c(rep(pca.out$rotation[11, 4], 2), rep(0, 2)), y = c(1, 0, 0, 1), col = "grey80")
abline(v = 0, lwd = 3)
axis(side = 2, at = c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), c("Plateletx1000perul", "RDW_Pct", "MCVfl", "PCVPct","Hgb.gperdL", "Fibrinogen", "EosinPct","MonoPct", "LymphPct","SegPct", "WBCx1000peruL"))
abline(v = .25, lty = 2, col = "grey20")
abline(v = -.25, lty = 2, col = "grey20")

# PCs by herd
#svg("./Plots/CBC_PCA_ByHerd_20160507.svg", height = 8, width = 10)
layout(matrix(c(1, 2, 2, 2, 3, 2, 2, 2), nrow = 2, byrow = T))
par(oma = c(1, 6, 1, 1), cex.lab = 1.2)
plot(1 ~ 1, cex = 0, xlim = c(-.6, .6), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC1", bty = "n", las = 2)
polygon(x = c(rep(pca.out$rotation[1, 1], 2), rep(0, 2)), y = c(11, 10, 10, 11), col = "grey80")
polygon(x = c(rep(pca.out$rotation[2, 1], 2), rep(0, 2)), y = c(10, 9, 9, 10), col = "grey80")
polygon(x = c(rep(pca.out$rotation[3, 1], 2), rep(0, 2)), y = c(9, 8, 8, 9), col = "grey80")
polygon(x = c(rep(pca.out$rotation[4, 1], 2), rep(0, 2)), y = c(8, 7, 7, 8), col = "grey80")
polygon(x = c(rep(pca.out$rotation[5, 1], 2), rep(0, 2)), y = c(7, 6, 6, 7), col = "grey80")
polygon(x = c(rep(pca.out$rotation[6, 1], 2), rep(0, 2)), y = c(6, 5, 5, 6), col = "grey80")
polygon(x = c(rep(pca.out$rotation[7, 1], 2), rep(0, 2)), y = c(5, 4, 4, 5), col = "grey80")
polygon(x = c(rep(pca.out$rotation[8, 1], 2), rep(0, 2)), y = c(4, 3, 3, 4), col = "grey80")
polygon(x = c(rep(pca.out$rotation[9, 1], 2), rep(0, 2)), y = c(3, 2, 2, 3), col = "grey80")
polygon(x = c(rep(pca.out$rotation[10, 1], 2), rep(0, 2)), y = c(2, 1, 1, 2), col = "grey80")
polygon(x = c(rep(pca.out$rotation[11, 1], 2), rep(0, 2)), y = c(1, 0, 0, 1), col = "grey80")
abline(v = 0, lwd = 3)
axis(side = 2, at = c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), las = 1,
     c("Plateletx1000perul", "RDW_Pct", "MCVfl", "PCVPct","Hgb.gperdL", "Fibrinogen", "EosinPct","MonoPct", "LymphPct","SegPct", "WBCx1000peruL"))
abline(v = .25, lty = 2, col = "grey20")
abline(v = -.25, lty = 2, col = "grey20")

plot(pca.out$x[, 2] ~ pca.out$x[, 1], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-5, 5), ylim = c(-5, 5),
     cex = cbc.small.complete$Movipneumonia.ELISA / 30, las = 1,
     xlab = "PC1", 
     ylab = "PC2")
leg.text <- levels(factor(cbc.small.complete$Herd))
legend("topleft", leg.text, pch = rep(16, 4), 
       col = seq(1:length(levels(factor(cbc.small.complete$Herd)))),
       ncol = 4, cex = 1.2)

plot(1 ~ 1, cex = 0, xlim = c(-.6, .6), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC2", bty = "n", las = 1)
polygon(x = c(rep(pca.out$rotation[1, 2], 2), rep(0, 2)), y = c(11, 10, 10, 11), col = "grey80")
polygon(x = c(rep(pca.out$rotation[2, 2], 2), rep(0, 2)), y = c(10, 9, 9, 10), col = "grey80")
polygon(x = c(rep(pca.out$rotation[3, 2], 2), rep(0, 2)), y = c(9, 8, 8, 9), col = "grey80")
polygon(x = c(rep(pca.out$rotation[4, 2], 2), rep(0, 2)), y = c(8, 7, 7, 8), col = "grey80")
polygon(x = c(rep(pca.out$rotation[5, 2], 2), rep(0, 2)), y = c(7, 6, 6, 7), col = "grey80")
polygon(x = c(rep(pca.out$rotation[6, 2], 2), rep(0, 2)), y = c(6, 5, 5, 6), col = "grey80")
polygon(x = c(rep(pca.out$rotation[7, 2], 2), rep(0, 2)), y = c(5, 4, 4, 5), col = "grey80")
polygon(x = c(rep(pca.out$rotation[8, 2], 2), rep(0, 2)), y = c(4, 3, 3, 4), col = "grey80")
polygon(x = c(rep(pca.out$rotation[9, 2], 2), rep(0, 2)), y = c(3, 2, 2, 3), col = "grey80")
polygon(x = c(rep(pca.out$rotation[10, 2], 2), rep(0, 2)), y = c(2, 1, 1, 2), col = "grey80")
polygon(x = c(rep(pca.out$rotation[11, 2], 2), rep(0, 2)), y = c(1, 0, 0, 1), col = "grey80")
abline(v = 0, lwd = 3)
axis(side = 2, at = c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), las = 1,
     c("Plateletx1000perul", "RDW_Pct", "MCVfl", "PCVPct","Hgb.gperdL", "Fibrinogen", "EosinPct","MonoPct", "LymphPct","SegPct", "WBCx1000peruL"))
abline(v = .25, lty = 2, col = "grey20")
abline(v = -.25, lty = 2, col = "grey20")
# dev.off()

# PCs by herd
par(mar = c(6, 6, 2, 2), oma = c(1, 1, 1, 1), cex.lab = 1, mfrow = c(1, 1))
plot(pca.out$x[, 3] ~ pca.out$x[, 4], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 30,
     xlab = "PC1: high for animals with high WBC, \n Fibrinogen, Hgb, PCV, MCV / low RDW", 
     ylab = "PC2: high for animals hight in LymphPct \n / low in SegPct, MonoPct")
leg.text <- levels(factor(cbc.small.complete$Herd))
legend("topleft", leg.text, pch = rep(16, 4), 
       col = seq(1:length(levels(factor(cbc.small.complete$Herd)))),
       ncol = 2)

# PCs by Movi qPCR
par(mar = c(6, 6, 2, 2), oma = c(1, 1, 1, 1), cex.lab = 1, mfrow = c(2, 3))
plot(pca.out$x[, 2] ~ pca.out$x[, 1], col = as.numeric(factor(cbc.small.complete$WADDL.Movi.PCR)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC1", 
     ylab = "PC2")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
# legend("topleft", leg.text, pch = rep(16, 4), 
#        col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
plot(pca.out$x[, 3] ~ pca.out$x[, 1], col = as.numeric(factor(cbc.small.complete$WADDL.Movi.PCR)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC1", 
     ylab = "PC3")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
# legend("topleft", leg.text, pch = rep(16, 4), 
#        col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
plot(pca.out$x[, 4] ~ pca.out$x[, 1], col = as.numeric(factor(cbc.small.complete$WADDL.Movi.PCR)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC1", 
     ylab = "PC4")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
# legend("topleft", leg.text, pch = rep(16, 4), 
#        col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
plot(pca.out$x[, 3] ~ pca.out$x[, 2], col = as.numeric(factor(cbc.small.complete$WADDL.Movi.PCR)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC2", 
     ylab = "PC3")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))

plot(pca.out$x[, 4] ~ pca.out$x[, 2], col = as.numeric(factor(cbc.small.complete$WADDL.Movi.PCR)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC2", 
     ylab = "PC4")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))

plot(pca.out$x[, 4] ~ pca.out$x[, 3], col = as.numeric(factor(cbc.small.complete$WADDL.Movi.PCR)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC3", 
     ylab = "PC4")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))

# PCs by Herd
par(mar = c(6, 6, 2, 2), oma = c(1, 1, 1, 1), cex.lab = 1, mfrow = c(2, 3))
plot(pca.out$x[, 2] ~ pca.out$x[, 1], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 8),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC1", 
     ylab = "PC2")
leg.text <- levels(factor(cbc.small.complete$Herd))
legend("topleft", leg.text, pch = rep(16, 4), ncol = 2, 
        col = seq(1:length(levels(factor(cbc.small.complete$Herd)))))
plot(pca.out$x[, 3] ~ pca.out$x[, 1], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC1", 
     ylab = "PC3")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
# legend("topleft", leg.text, pch = rep(16, 4), 
#        col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
plot(pca.out$x[, 4] ~ pca.out$x[, 1], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC1", 
     ylab = "PC4")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
# legend("topleft", leg.text, pch = rep(16, 4), 
#        col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
plot(pca.out$x[, 3] ~ pca.out$x[, 2], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC2", 
     ylab = "PC3")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))

plot(pca.out$x[, 4] ~ pca.out$x[, 2], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC2", 
     ylab = "PC4")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))

plot(pca.out$x[, 4] ~ pca.out$x[, 3], col = as.numeric(factor(cbc.small.complete$Herd)), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6),
     cex = cbc.small.complete$Movipneumonia.ELISA / 40,
     xlab = "PC3", 
     ylab = "PC4")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))

#---------------------------#
#---- quick LDA test -------#
#---------------------------#
require(MASS)
lda.fit1 <- lda(qPCRResult ~ WBCx1000peruL + SegPct + LymphPct + MonoPct +
                  EosinPct + Fibrinogen + Hgb.gperdL + PCVPct + MCVfl + RDW_Pct + Plateletx1000peruL,
                data = cbc.small.complete)
# plot LDA loadings
par(mfrow = c(1, 1), las = 1, cex.axis = .8, oma = c(0, 2, 0, 0), mar = c(6, 6, 1, 1))
plot(1 ~ 1, cex = 0, xlim = c(-.2, .2), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC1", bty = "n")
polygon(x = c(rep(lda.fit1$scaling[1, 1], 2), rep(0, 2)), y = c(11, 10, 10, 11), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[2, 1], 2), rep(0, 2)), y = c(10, 9, 9, 10), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[3, 1], 2), rep(0, 2)), y = c(9, 8, 8, 9), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[4, 1], 2), rep(0, 2)), y = c(8, 7, 7, 8), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[5, 1], 2), rep(0, 2)), y = c(7, 6, 6, 7), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[6, 1], 2), rep(0, 2)), y = c(6, 5, 5, 6), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[7, 1], 2), rep(0, 2)), y = c(5, 4, 4, 5), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[8, 1], 2), rep(0, 2)), y = c(4, 3, 3, 4), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[9, 1], 2), rep(0, 2)), y = c(3, 2, 2, 3), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[10, 1], 2), rep(0, 2)), y = c(2, 1, 1, 2), col = "grey80")
polygon(x = c(rep(lda.fit1$scaling[11, 1], 2), rep(0, 2)), y = c(1, 0, 0, 1), col = "grey80")
abline(v = 0, lwd = 3)
axis(side = 2, at = c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5), 
     c("Plateletx1000perul", "RDW_Pct", "MCVfl", "PCVPct","Hgb.gperdL", "Fibrinogen", "EosinPct","MonoPct", "LymphPct","SegPct", "WBCx1000peruL"))
abline(v = .25, lty = 2, col = "grey20")
abline(v = -.25, lty = 2, col = "grey20")

lda.preds <- predict(lda.fit1)
fit.check <- cbind(lda.preds$class, as.character(cbc.small.complete$WADDL.Movi.PCR))
table(fit.check[, 1], fit.check[ ,2])

# PCs by Movi qPCR, split by age
complete.juvs.f <- subset(cbc.small.complete, Sex == "F" & (ACLASS_SAMP == "Lamb" | ACLASS_SAMP == "Yearling"))
complete.juvs.m <- subset(cbc.small.complete, Sex == "M" & (ACLASS_SAMP == "Lamb" | ACLASS_SAMP == "Yearling"))
complete.ads.f <- subset(cbc.small.complete, Sex == "F" & (ACLASS_SAMP == "Adult"))
complete.ads.m <- subset(cbc.small.complete, Sex == "M" & (ACLASS_SAMP == "Adult"))
par(mar = c(6, 6, 2, 2), oma = c(1, 1, 1, 1), cex.lab = 1, mfrow = c(2, 2))
plot(complete.juvs.f$PC2 ~ complete.juvs.f$PC1, col = as.numeric(factor(complete.juvs.f$WADDL.Movi.PCR, levels = c("Detected", "Indeterminate", "Not detected"))), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6), main = "Juvenile females",
     cex = complete.juvs.f$Movipneumonia.ELISA / 30,
     xlab = "PC1: high for animals with high WBC, \n Fibrinogen, Hgb, PCV, MCV / low RDW", 
     ylab = "PC2: high for animals hight in LymphPct \n / low in SegPct, MonoPct")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
legend("topleft", leg.text, pch = rep(16, 4), 
       col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
# young rams
plot(complete.juvs.m$PC2 ~ complete.juvs.m$PC1, 
     col = as.numeric(factor(complete.juvs.m$WADDL.Movi.PCR, levels = c("Detected", "Indeterminate", "Not detected"))), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6), main = "Juvenile males",
     cex = complete.juvs.m$Movipneumonia.ELISA / 30,
     xlab = "PC1: high for animals with high WBC, \n Fibrinogen, Hgb, PCV, MCV / low RDW", 
     ylab = "PC2: high for animals hight in LymphPct \n / low in SegPct, MonoPct")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
legend("topleft", leg.text, pch = rep(16, 4), 
       col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
# adult ewes
plot(complete.ads.f$PC2 ~ complete.ads.f$PC1, col = as.numeric(factor(complete.ads.f$WADDL.Movi.PCR, levels = c("Detected", "Indeterminate", "Not detected"))), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6), main = "Adult females",
     cex = complete.ads.f$Movipneumonia.ELISA / 30,
     xlab = "PC1: high for animals with high WBC, \n Fibrinogen, Hgb, PCV, MCV / low RDW", 
     ylab = "PC2: high for animals hight in LymphPct \n / low in SegPct, MonoPct")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
legend("topleft", leg.text, pch = rep(16, 4), 
       col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))))
# adult rams
plot(complete.ads.m$PC2 ~ complete.ads.m$PC1, 
     col = as.numeric(factor(complete.ads.m$WADDL.Movi.PCR, levels = c("Detected", "Indeterminate", "Not detected"))), 
     pch = 16, xlim = c(-6, 6), ylim = c(-6, 6), main = "Adult males",
     cex = complete.ads.m$Movipneumonia.ELISA / 30,
     xlab = "PC1: high for animals with high WBC, \n Fibrinogen, Hgb, PCV, MCV / low RDW", 
     ylab = "PC2: high for animals hight in LymphPct \n / low in SegPct, MonoPct")
leg.text <- levels(factor(cbc.small.complete$WADDL.Movi.PCR))
legend("topleft", leg.text, pch = rep(16, 4), 
       col = seq(1:length(levels(factor(cbc.small.complete$WADDL.Movi.PCR)))), ncol = 3)

# by herd
aso.cbc <- subset(cbc.small.complete, Herd == "Asotin")
los.cbc <- subset(cbc.small.complete, Herd == "Lostine")
imn.cbc <- subset(cbc.small.complete, Herd == "Imnaha")
bc.cbc <- subset(cbc.small.complete, Herd == "Big Canyon")

par(mfrow = c(2, 2))
plot(aso.cbc$PC2 ~ aso.cbc$PC1, col = aso.cbc$qPCRResult, pch = 16,
     cex = aso.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))
plot(bc.cbc$PC2 ~ bc.cbc$PC1, col = bc.cbc$qPCRResult, pch = 16,
     cex = bc.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))
plot(imn.cbc$PC2 ~ imn.cbc$PC1, col = imn.cbc$qPCRResult, pch = 16,
     cex = imn.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))
plot(los.cbc$PC2 ~ los.cbc$PC1, col = los.cbc$qPCRResult, pch = 16,
     cex = los.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))

# colored by age class
par(mfrow = c(2, 2))
plot(aso.cbc$PC2 ~ aso.cbc$PC1, col = factor(aso.cbc$ACLASS_SAMP, levels = levels(cbc.small.complete$ACLASS_SAMP)), pch = 16,
     cex = aso.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))
plot(bc.cbc$PC2 ~ bc.cbc$PC1, col = factor(bc.cbc$ACLASS_SAMP, levels = levels(cbc.small.complete$ACLASS_SAMP)), pch = 16,
     cex = bc.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))
plot(imn.cbc$PC2 ~ imn.cbc$PC1, col = factor(imn.cbc$ACLASS_SAMP, levels = levels(cbc.small.complete$ACLASS_SAMP)), pch = 16,
     cex = imn.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))
plot(los.cbc$PC2 ~ los.cbc$PC1, col = factor(los.cbc$ACLASS_SAMP, levels = levels(cbc.small.complete$ACLASS_SAMP)), pch = 16,
     cex = los.cbc$Movipneumonia.ELISA / 25, xlim = c(-4, 4), ylim = c(-4, 4))

# CBC by PCR status
par(mfrow = c(3, 4), oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 2))
# plot(as.numeric(as.character(cbc.data$Movipneumonia.ELISA)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
#      type = "p", col = "black", xaxt = "n", xlab = "Movi PCR", ylab = "Movi ELISA (% Inhibition)", las = 1, xlim = c(0.5, 3.5))
# abline(h = 40, lty = 2)
plot(as.numeric(as.character(cbc.data$WBCx1000peruL)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", ylab = "WBC", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$SegPct)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", ylab = "Seg %", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))

plot(as.numeric(as.character(cbc.data$LymphPct)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "Lymph %", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))

plot(as.numeric(as.character(cbc.data$MonoPct)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "Mono %", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))

plot(as.numeric(as.character(cbc.data$Fibrinogen)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "Fibrinogen", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$EosinPct)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "Eosinophil %", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$Fibrinogen)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "Fibrinogen", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$Hgb.gperdL)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "Hemoglobin (g/dL)", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$PCVPct)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "PVC", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$MCV)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "MCV", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$RDW)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "RDW", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))
plot(as.numeric(as.character(cbc.data$Plateletx1000peruL)) ~ jitter(as.numeric(factor(cbc.data$WADDL.Movi.PCR)), .2), 
     type = "p", col = rgb(0, 0, 0, alpha = .5), xaxt = "n", xlab = "Movi PCR", 
     ylab = "Platelets per 1000 uL", las = 1, xlim = c(0.5, 3.5))
axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))

#abline(h = 40, lty = 2)


#axis(side = 1, labels = levels(factor(full.data$WADDL.Movi.PCR)), at = c(1:3))


