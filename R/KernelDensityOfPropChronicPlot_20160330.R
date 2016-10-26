full.data <- read.csv("./Data/LostineMoviPastLung_150630.csv", header = T)
full.data <- read.csv("./Data/LostineSampleData_20160406_v2.csv", header = T)
full.data$id <- full.data$Animal.ID
full.data$movi_qpcr <- full.data$WADDL.Movi.PCR
full.data$cap_bioyr <- full.data$CAP_BIOYR
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")


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

# reduced to just animals with 3 or more samples
data <- subset(full.data , !(id %in% observed.once | id %in% observed.twice))
data$id <- factor(data$id)

pcr.by.id.full <- table(full.data$qPCRResult, full.data$id)
prop.pos.full <- prop.table(pcr.by.id.full, margin = 2)[1, ]

pcr.by.id <- table(data$qPCRResult, data$id)
prop.pos <- prop.table(pcr.by.id, margin = 2)[1, ]
rbind(pcr.by.id, 1-prop.pos)

par(mfrow = c(1, 1), mar = c(4, 4, 2, 4))
hist(1 - prop.pos, breaks = 10, col = rgb(0, 0, 0, alpha = .5),
     xlim = c(0, 1), ylim = c(0, 10),
     main = "", ylab = "# Animals", 
     xlab = "Proportion positive", 
     las = 1)
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
mtext(side = 4, line = 3, "Density")
axis(side = 4, at = seq(0, max(d$y), by = .25), labels = seq(0, max(d$y), by = .25), las = 1)
leg.text <- "Density over all tested animals"
legend(y = max(d$y) + .2, x = -.05, leg.text, lty = 1, bty = "n")
