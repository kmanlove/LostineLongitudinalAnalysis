compd.data <- read.csv("./Data/compiled_data_summary_151115b.csv", header = T)
los <- subset(compd.data, Pop == "Lostine")
los$CLASS[43:44] <- "LAMBS"

samples <- read.csv("./Data/LostineSampleData_20160406.csv", header = T)

# strip off sampled animals
sampled.anims <- levels(factor(samples$Animal.ID))
cohort <- rep(NA, length(sampled.anims))
for(i in 1:length(sampled.anims)){
  k <- subset(samples, as.character(Animal.ID) == as.character(levels(factor(samples$Animal.ID))[i]))
  cohort[i] <- k$BirthYear[1]
}

levels(factor(samples$Animal.ID))[cohort == 2010]

cohort.fac <- factor(cohort, levels = c(seq(1990, 2014)))
sampled.anims.by.cohort <- as.vector(table(cohort.fac))
los.short <- subset(los, year %in% c(1990:2014))

plot(sampled.anims.by.cohort ~ los.short$Lambs, pch = 16,
     xlim = c(0, 25), ylim = c(0, 14),
     xlab = "Number recruits",
     ylab = "Number sampled animals in cohort")
text(sampled.anims.by.cohort ~ los.short$Lambs, 
     labels = as.character(seq(1990,2014)),
     pos = 3, cex = .8)
lm.test <- lm(sampled.anims.by.cohort ~ los.short$Lambs)
abline(a = 0, b = 1, lty = 2)
