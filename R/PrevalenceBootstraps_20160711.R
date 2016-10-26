# Movi prevalence estimation
samples <- read.csv("./Data/LostineSample2016Final_WithClasses_20160927.csv", header = T)
samples <- subset(samples, movi_qpcr %in% c("Detected", "Not detected"))
samples$qPCRResult <- ifelse(samples$movi_qpcr == "Detected" | samples$movi_qpcr == "Indeterminate", "P", "N")
samples <- subset(samples, cap_bioyr >= 2008)

#samples <- read.csv("./Data/LostineSample2016Final.csv", header = T)
#length(levels(factor(samples$id)))

# prev based on first sampling event for EVERYBODY
# strip out each animal's first sampling event
first.event <- vector("list", length(levels(factor(samples$id))))
for(i in 1:length(levels(factor(samples$id)))){
  k <- subset(samples, as.character(id) == as.character(levels(factor(samples$id))[i]) & cap_bioyr >= 2011)
  first.event[[i]] <- k[which.min(as.numeric(as.Date(k$capture_date, format = "%m/%d/%Y"))), ]
}

first.event.frame <- do.call("rbind", first.event)
hist(first.event.frame$cap_bioyr, col = "grey60")
length(levels(factor(first.event.frame$id)))
table(first.event.frame$movi_qpcr)

# subset adults/yrs/lambs
first.event.ads <- subset(first.event.frame, age_class == "Adult" & sex == "F")
first.event.lamb <- subset(first.event.frame, age_class == "Lamb")
first.event.yrl <- subset(first.event.frame, age_class == "Yearling")
first.event.frame$juv.ind <- ifelse(first.event.frame$age_class == "Adult", "adult", "juv")
table(first.event.ads$movi_qpcr)
table(first.event.juv$movi_qpcr)
 
# bootstrapped prevalence across all years in adult females
nboot <- 1000
boot.ewe.prev <- rep(NA, nboot)
boot.lamb.prev <- rep(NA, nboot)
boot.yrl.prev <- rep(NA, nboot)
for(i in 1:nboot){
  k <- first.event.ads[sample(1:dim(first.event.ads)[1], size = dim(first.event.ads)[1], replace = T), ]
  boot.ewe.prev[i] <- table(k$movi_qpcr, useNA = "always")[1]/dim(first.event.ads)[1]
}

for(i in 1:nboot){
  k <- first.event.lamb[sample(1:dim(first.event.lamb)[1], size = dim(first.event.lamb)[1], replace = T), ]
  boot.lamb.prev[i] <- table(k$movi_qpcr, useNA = "always")[1]/dim(first.event.lamb)[1]
}

for(i in 1:nboot){
  k <- first.event.yrl[sample(1:dim(first.event.yrl)[1], size = dim(first.event.yrl)[1], replace = T), ]
  boot.yrl.prev[i] <- table(k$movi_qpcr, useNA = "always")[1]/dim(first.event.yrl)[1]
}

adewe.prev.all.yrs <- quantile(boot.ewe.prev, c(0.025, 0.975))
lamb.prev.all.yrs <- quantile(boot.lamb.prev, c(0.025, 0.975))
yrl.prev.all.yrs <- quantile(boot.yrl.prev, c(0.025, 0.975))

# Comparison of WITHIN-YEAR prevalence
# need first sampling event of each animal in each year

prev.xtab <- table(first.event.frame$age_class, first.event.frame$movi_qpcr)
age.spec.prevs <- prev.xtab[, 1] / table(first.event.frame$age_class)

require(graphics)
#svg("./Plots/AgeSpecPrev_20160711.svg", height = 4, width = 3)
plot(x = 0, y = 0, xlim = c(.5, 3.5), ylim = c(0, 1), cex = 0,
     xaxt = "n", ylab = "Prevalence", las = 1, xlab = "")
segments(x0 = 1, x1 = 1, y0 = adewe.prev.all.yrs[1], y1 = adewe.prev.all.yrs[2], lwd = 2)
segments(x0 = 2, x1 = 2, y0 = yrl.prev.all.yrs[1], y1 = yrl.prev.all.yrs[2], lwd = 2)
segments(x0 = 3, x1 = 3, y0 = lamb.prev.all.yrs[1], y1 = lamb.prev.all.yrs[2], lwd = 2)
segments(x0 = 0.9, x1 = 1.1, y0 = age.spec.prevs[1], y1 = age.spec.prevs[1], lwd = 2)
segments(x0 = 1.9, x1 = 2.1, y0 = age.spec.prevs[3], y1 = age.spec.prevs[3], lwd = 2)
segments(x0 = 2.9, x1 = 3.1, y0 = age.spec.prevs[2], y1 = age.spec.prevs[2], lwd = 2)
axis(side = 1, at = c(1, 2, 3), labels = c("Adult Ewes", "Yearlings", "Lambs"))
#dev.off()
