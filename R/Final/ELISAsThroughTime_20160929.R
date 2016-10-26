full.data <- read.csv("./Data/LostineSampleData_20160505_v1.csv", header = T)
full.data$id <- full.data$Animal.ID
full.data$movi_qpcr <- full.data$WADDL.Movi.PCR
full.data$cap_bioyr <- full.data$CAP_BIOYR
full.data <- subset(full.data, movi_qpcr != "Indeterminate")
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")
full.data$CaptureDateContinuous <- difftime((as.POSIXlt(strptime(full.data$Capture.Date, format = "%d/%m/%Y"))), 
                                            min(as.POSIXlt(strptime(full.data$Capture.Date, format = "%d/%m/%Y"))),
                                            units = "days")

boxplot(full.data$Movi.ELISA ~ full.data$Age.Class)
table(full.data$id, full.data$Age.Class)
plot(full.data$Movi.ELISA ~ full.data$CaptureDateContinuous, col = full.data$Age.Class)

full.data.recent <- subset(full.data, CaptureDateContinuous >= 4000)
ELISA.fit <- lm(full.data.recent$Movi.ELISA ~ full.data.recent$CaptureDateContinuous + full.data.recent$Age.Class)
summary(ELISA.fit)

plot(full.data$Movi.ELISA ~ full.data$CaptureDateContinuous, ylab = "Movi ELISA",
     xlab = "Sampling date",
     col = full.data$Age.Class, xlim = c(4500, 6000), xaxt = "n")
axis(side = 1, at = c(4800, 4800 + 365, 4800 + 365*2, 4800+365*3),
     labels = c("2012-2013", "2013-2014", "2014-2015", "2015-2016"))

kruskal.test(full.data$Movi.ELISA ~ factor(full.data$qPCRResult))
