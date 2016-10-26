#---------------------------------------------#
#-- BHS vs. Domestic qPCR and ELISA data -----#
#---------------------------------------------#


#-- Read in and prep dataframe --#

#elisa.data <- read.csv(file.choose(), header = T, sep = "\t")
elisa.data <- read.csv("~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Data/BHSLongitudinal_ELISA_qPCR_23Sept2014.csv", header = T)
#-- csv is long.sampling.data.csv in work/kez/research/ecolpapers/lostinelongsampling/data/long.sampling.data.csv ish --#
#-- read in the csv version of the longitudinal sampling data. The file.choose part of the command below opens a navigation window that --#
#-- will let you navigate directly to the dataset --#

head(elisa.data) 

qpcr.data <- read.csv("~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Data/BHS_qPCRData_23Sept2014.csv", header = T)

elisa.data$qPCR <- rep(NA, dim(elisa.data)[1])
for(i in 1:dim(elisa.data)[1]){
  k <- subset(qpcr.data, as.character(Animal_ID) == as.character(elisa.data$Animal_ID[i]) & as.character(CaptureDate) == as.character(elisa.data$CaptureDate[i]))
  elisa.data$qPCR[i] <- ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$CT_Magmax[1])))
}

full.bhs.dat <- subset(elisa.data, is.na(qPCR) == F & is.na(Movipneumonia.ELISA) == F)
plot(full.bhs.dat$Movipneumonia.ELISA ~ full.bhs.dat$qPCR, xlim = c(0, 40), ylim = c(0, 100), pch = 16, ylab = "Movi ELISA", xlab = "Movi qPCR (CT_Magmax)")

#-- Shows the first five rows of the dataset. I use this to be sure that the dataset was read in with the correct delimiting.
elisa.data$WADDL.Movi.PCR <- ifelse(elisa.data$WADDL.Movi.PCR == "Indeterminant", "Indeterminate", as.character(elisa.data$WADDL.Movi.PCR))
elisa.data$WADDL.Movi.PCR <- factor(elisa.data$WADDL.Movi.PCR)
#-- These two lines correct the indeterminate/indeterminant inconsistency in the dataset. 
#-- After these lines, all indeterminate shedders are labeled "Indeterminate".
table(elisa.data$WADDL.Movi.PCR)

#-- build and explore determination cohort --#

Animal.tab <- table(elisa.data$Animal_ID)
#-- Tables number of observations for each Animal ID and stores table in "Animal.tab" object
AnimalByOutcome.tab <- table(elisa.data$Animal_ID,elisa. data$WADDL.Movi.PCR)
Animal.tab[which(Animal.tab >= 4)]
#-- Extracts Animal IDs with 4 or more samples
Animal.tab[which(Animal.tab == 3)]
#-- Extracts Animal IDs with 3 samples

ManySamplesOutcomes.tab <- AnimalByOutcome.tab[which(Animal.tab >= 4), ]
ManySamplesProps <- ManySamplesOutcomes.tab[ ,1] / (ManySamplesOutcomes.tab[ ,1] + ManySamplesOutcomes.tab[ ,4])

#--------------------------------------------------------------------------------#
#-- Histogram of proportion of positive tests for sheep with 4 or more samples --#
#--------------------------------------------------------------------------------#

#hist(ManySamplesProps, breaks = 10, col = "grey80", xlab = "Proportion of PCR-positive tests for \n the 12 animals tested 4 or more times", main = "")
#-- build a histogram of the proportion of tests that were positive for animals tested 4 or more times

#-- classify animals in determination cohort into carriage states ----#
#-- place cut line at .6: anybody above .6 is classified as chronic --#
determination.cohort <- subset(elisa.data, Animal_ID %in% names(which(Animal.tab >= 4)))
table(factor(determination.cohort$Animal_ID), determination.cohort$WADDL.Movi.PCR)
Chronic <- c("01LO34", "04LO12", "99L09")
NearlyAlwaysHealthy <- c("04LO60", "04LO73", "05LO68", "11LO40", "11LO44", "12LO56", "10LO33", "08LO29")

full.bhs.dat$AnimalCarriage <- ifelse(full.bhs.dat$Animal_ID %in% Chronic, "red", ifelse(full.bhs.dat$Animal_ID %in% NearlyAlwaysHealthy, "blue", "grey70"))

plot(full.bhs.dat$Movipneumonia.ELISA ~ log(40 - full.bhs.dat$qPCR), xlim = c(0, 3.2), ylim = c(0, 100), pch = 16, col = full.bhs.dat$AnimalCarriage, cex = 2)
plot(full.bhs.dat$Movipneumonia.ELISA ~ full.bhs.dat$qPCR, xlim = c(0, 40), ylim = c(0, 100), pch = 16, col = full.bhs.dat$AnimalCarriage, cex = 2)



determination.cohort$CarriageState <- ifelse(determination.cohort$Animal_ID %in% Chronic, "chronic", "not.chronic")

PropPos <- matrix(NA, ncol = 4, nrow = length(levels(factor(determination.cohort$Animal_ID))))
TrueStatus <- rep(NA, length(levels(factor(determination.cohort$Animal_ID))))
for(i in 1:4){
  tests.to.use <- subset(determination.cohort, TestNumber <= i)
  for(j in 1:length(levels(factor(determination.cohort$Animal_ID)))){
    anim.spec.data <- subset(tests.to.use, Animal_ID == levels(factor(determination.cohort$Animal_ID))[j])
    anim.spec.results <- table(anim.spec.data$WADDL.Movi.PCR)
    PropPos[j, i] <- anim.spec.results[1] / i
    TrueStatus[j] <- anim.spec.data$CarriageState[1]
  }
}
