data <- read.csv("./Data/TestResultsOverYears_DataToModel_20160926_R.csv", header = T)

all.pos.yr1 <- subset(data, AllPos == 1)
table(all.pos.yr1$AnyPos.1)
allpos1.anypos2 <- prop.test(x = table(all.pos.yr1$AnyPos.1))
1 - allpos1.anypos2$conf.int

no.pos.yr1 <- subset(data, AnyPos == 0)
nopos1.nopos2 <- prop.test(x = table(no.pos.yr1$AnyPos.1))
nopos1.nopos2$conf.int
