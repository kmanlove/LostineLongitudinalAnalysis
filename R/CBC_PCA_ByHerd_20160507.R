cbc.data <- read.csv("./Data/BHS_CBC_Data_V2.csv", sep = ",")


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


# PCs by herd
#svg("./Plots/CBC_PCA_ByHerd_20160509_v2.svg", height = 8, width = 10)
#pt.cols <- rainbow(n = 13)

pt.cols <- c("cyan2", "green", "darkviolet", "navyblue", "deeppink1", "blue", "grey60", "brown", "forestgreen", 
             "gold",  "red", "lightblue", "black")
#pt.cols <- cm.colors(n = 13)
layout(matrix(c(1, 2, 2, 2, 3, 2, 2, 2), nrow = 2, byrow = T))
par(oma = c(1, 6, 1, 1), cex.lab = 1.2)
plot(1 ~ 1, cex = 0, xlim = c(-.6, .6), ylim = c(0, 11.5), ylab = "", yaxt = "n", 
     xlab = "Correlation with PC1", bty = "n", las = 1)
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

plot(pca.out$x[, 2] ~ pca.out$x[, 1], col = pt.cols[as.numeric(factor(cbc.small.complete$Herd))], 
     pch = 16, xlim = c(-5, 5), ylim = c(-5, 5),
#     cex = cbc.small.complete$Movipneumonia.ELISA / 30, 
     cex = 2.5, 
     las = 1,
     xlab = "PC1", 
     ylab = "PC2", cex.lab = 1.2, cex.axis = 1.2)
leg.text <- levels(factor(cbc.small.complete$Herd))
legend("top", leg.text, pch = rep(16, 4), 
       col = pt.cols[seq(1:length(levels(factor(cbc.small.complete$Herd))))],
       ncol = 4, pt.cex = 2, bty = "n", cex = 1.4)
abline(h = 4.0)

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
#dev.off()