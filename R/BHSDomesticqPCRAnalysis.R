# load the lme4 package
require(lme4)

# read in Bighorn_Sheep_QPCR_141026 dataset
bhs.data <- read.csv("~/work/Kezia/Research/EcologyPapers/LostineLongitudinalStudy/Data/Bighorn_Sheep_QPCR_141026.csv", header = T)
names(bhs.data)

# reset dates to be R date class
bhs.data$Sample1Date <- as.character(as.vector(strptime(bhs.data$Sample1, format = "%m/%d/%Y")))
bhs.data$Sample2Date <- as.character(as.vector(strptime(bhs.data$Sample2, format = "%m/%d/%Y")))
bhs.data$Sample3Date <- as.character(as.vector(strptime(bhs.data$Sample3, format = "%m/%d/%Y")))
bhs.data$Sample4Date <- as.character(as.vector(strptime(bhs.data$Sample4, format = "%m/%d/%Y")))
bhs.data$Sample5Date <- as.character(as.vector(strptime(bhs.data$Sample5, format = "%m/%d/%Y")))

bhs.data$CtDirectQIAmp1.Detected <- ifelse(bhs.data$CT_Direct_NW_QIAmp == "ND", 0, ifelse(is.na(bhs.data$CT_Direct_NW_QIAmp) == F, 1, NA))
bhs.data$CtDirectQIAmp2.Detected <- ifelse(bhs.data$CT_Direct_NW_QIAmp2 == "ND", 0, ifelse(is.na(bhs.data$CT_Direct_NW_QIAmp2) == F, 1, NA))
bhs.data$CtDirectQIAmp3.Detected <- ifelse(bhs.data$CT_Direct_NW_QIAmp3 == "ND", 0, ifelse(is.na(bhs.data$CT_Direct_NW_QIAmp3) == F, 1, NA))
bhs.data$CtDirectQIAmp4.Detected <- ifelse(bhs.data$CT_Direct_NW_QIAmp4 == "ND", 0, ifelse(is.na(bhs.data$CT_Direct_NW_QIAmp4) == F, 1, NA))
bhs.data$CtDirectQIAmp5.Detected <- ifelse(bhs.data$CT_Direct_NW_QIAmp5 == "ND", 0, ifelse(is.na(bhs.data$CT_Direct_NW_QIAmp5) == F, 1, NA))

data.check <- subset(bhs.data, select = c("CT_Direct_NW_QIAmp", "CtDirectQIAmp1.Detected"))

# unpack rows into columns
StartDate <- StopDate <- AnimalID <- ResultIn <- ResultOut <- CtDirectSwabIn <- CtDirectSwabOut <- CtDirectQIAmpIn <- CtDirectQIAmpOut <- CtDirectQIAmpBinaryIn <- CtDirectQIAmpBinaryOut <- rep(NA, 4 * dim(bhs.data)[1])
for(i in 1:dim(bhs.data)[1]){
  StartDate[i] <- bhs.data$Sample1Date[i]
  StartDate[(2 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Sample2Date[i]
  StartDate[(3 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Sample3Date[i]
  StartDate[(4 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Sample4Date[i]

  StopDate[i] <- bhs.data$Sample2Date[i]
  StopDate[(2 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Sample3Date[i]
  StopDate[(3 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Sample4Date[i]
  StopDate[(4 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Sample5Date[i]

  AnimalID[i] <- bhs.data$AnimalID[i]
  AnimalID[(2 - 1) * dim(bhs.data)[1] + i] <- bhs.data$AnimalID[i]
  AnimalID[(3 - 1) * dim(bhs.data)[1] + i] <- bhs.data$AnimalID[i]
  AnimalID[(4 - 1) * dim(bhs.data)[1] + i] <- bhs.data$AnimalID[i]

  ResultIn[i] <- bhs.data$Result1[i]
  ResultIn[(2 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Result2[i]
  ResultIn[(3 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Result3[i]
  ResultIn[(4 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Result4[i]

  ResultOut[i] <- bhs.data$Result2[i]
  ResultOut[(2 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Result3[i]
  ResultOut[(3 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Result4[i]
  ResultOut[(4 - 1) * dim(bhs.data)[1] + i] <- bhs.data$Result5[i]

  CtDirectSwabIn[i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM[i]))
  CtDirectSwabIn[(2 - 1) * dim(bhs.data)[1] + i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM2[i]))
  CtDirectSwabIn[(3 - 1) * dim(bhs.data)[1] + i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM3[i]))
  CtDirectSwabIn[(4 - 1) * dim(bhs.data)[1] + i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM4[i]))
  
  CtDirectSwabOut[i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM2[i]))
  CtDirectSwabOut[(2 - 1) * dim(bhs.data)[1] + i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM3[i]))
  CtDirectSwabOut[(3 - 1) * dim(bhs.data)[1] + i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM4[i]))
  CtDirectSwabOut[(4 - 1) * dim(bhs.data)[1] + i] <- as.numeric(as.character(bhs.data$Ct_Direct_swab_MM5[i]))

  CtDirectQIAmpIn[i] <- as.numeric(as.character(bhs.data$CT_Direct_NW_QIAmp[i]))
  CtDirectQIAmpIn[(2 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CT_Direct_NW_QIAmp2[i])
  CtDirectQIAmpIn[(3 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CT_Direct_NW_QIAmp3[i])
  CtDirectQIAmpIn[(4 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CT_Direct_NW_QIAmp4[i])
  
  CtDirectQIAmpOut[i] <- as.numeric(as.character(bhs.data$CT_Direct_NW_QIAmp2[i]))
  CtDirectQIAmpOut[(2 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CT_Direct_NW_QIAmp3[i])
  CtDirectQIAmpOut[(3 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CT_Direct_NW_QIAmp4[i])
  CtDirectQIAmpOut[(4 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CT_Direct_NW_QIAmp5[i])

  CtDirectQIAmpBinaryIn[i] <- as.numeric(as.character(bhs.data$CtDirectQIAmp1.Detected[i]))
  CtDirectQIAmpBinaryIn[(2 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CtDirectQIAmp2.Detected[i])
  CtDirectQIAmpBinaryIn[(3 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CtDirectQIAmp3.Detected[i])
  CtDirectQIAmpBinaryIn[(4 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CtDirectQIAmp4.Detected[i])
  
  CtDirectQIAmpBinaryOut[i] <- as.numeric(as.character(bhs.data$CtDirectQIAmp2.Detected[i]))
  CtDirectQIAmpBinaryOut[(2 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CtDirectQIAmp3.Detected[i])
  CtDirectQIAmpBinaryOut[(3 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CtDirectQIAmp4.Detected[i])
  CtDirectQIAmpBinaryOut[(4 - 1) * dim(bhs.data)[1] + i] <- as.character(bhs.data$CtDirectQIAmp5.Detected[i])
}

bhs.long <- as.data.frame(cbind(StartDate, StopDate, ResultIn, ResultOut, CtDirectSwabIn, CtDirectSwabOut, CtDirectQIAmpIn, CtDirectQIAmpOut))
bhs.long$Transition <- paste(bhs.long$ResultIn, bhs.long$ResultOut, sep = "")
bhs.long.nonas <- subset(bhs.long, is.na(ResultIn) == F & is.na(ResultOut) == F)
table(bhs.long.nonas$Transition)
bhs.wide.check <- subset(bhs.data, select = c("AnimalID", "Result1", "Result2", "Result3", "Result4", "Result5", "Result6"))
StartDateDs <- StopDateDs <- ResultInDs <- ResultOutDs <- rep(NA, dim(ds.qpcr)[1] * 2)

data.PosTrial1 <- bhs.long
data.PosTrial1$QIAmp.NoND <- ifelse(data.PosTrial1$CtDirectQIAmpOut == "ND", 40, as.numeric(as.character(data.PosTrial1$CtDirectQIAmpOut)))
data.PosTrial1$QIAmp.NoND.In <- ifelse(data.PosTrial1$CtDirectQIAmpIn == "ND", 40, as.numeric(as.character(data.PosTrial1$CtDirectQIAmpIn)))
data.PosTrial1$Result2.binary <- (ifelse(data.PosTrial1$ResultOut == 2, NA, data.PosTrial1$ResultOut)) - 1

#-----------------------------------------------------------------------------------#
#-- Sample sizes -------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
tot.enriched <- 
tot.enriched.pairs <- table(!(is.na(data.PosTrial1$ResultIn)), !(is.na(data.PosTrial1$ResultOut)))
enriched.pairs <- subset(data.PosTrial1, !(is.na(ResultIn)) & !(is.na(ResultOut)))
table(enriched.pairs$ResultIn, enriched.pairs$ResultOut)

tot.QIAmp.pairs <- table(!(is.na(data.PosTrial1$QIAmp.NoND.In)), !(is.na(data.PosTrial1$QIAmp.NoND)))
QIAmp.pairs <- subset(data.PosTrial1, !(is.na(QIAmp.NoND.In)) & !(is.na(QIAmp.NoND)))

cbind(data.PosTrial1$QIAmp.NoND.In, data.PosTrial1$QIAmp.NoND)

#-----------------------------------------------------------------------------------#
#-- Is qPCR result (either swab or wash) at time t predictive of result at t + 1? --#
#-----------------------------------------------------------------------------------#
result.preds.result <- glm(Result2.binary ~ ResultIn, family = binomial(link = "logit"), data = data.PosTrial1)
QIAmp.preds.result <- glm(Result2.binary ~ as.numeric(as.character(QIAmp.NoND)), family = binomial(link = "logit"), data = data.PosTrial1)
QIAmp.preds.QIAmp <- lm(as.numeric(as.character(QIAmp.NoND)) ~ as.numeric(as.character(QIAmp.NoND.In)), data = data.PosTrial1)
swab.preds.result <- glm(Result2.binary ~ as.numeric(as.character(CtDirectSwabIn)), family = binomial(link = "logit"), data = data.PosTrial1)
swab.preds.QIAmp <- lm(as.numeric(as.character(QIAmp.NoND)) ~ as.numeric(as.character(CtDirectSwabIn)), data = data.PosTrial1)

summary(QIAmp.preds.result)
summary(swab.preds.result)
summary(result.preds.result)

