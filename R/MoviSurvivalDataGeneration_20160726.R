full.data <- read.csv("./Data/LostineSampleData_20160505_v1.csv", header = T)
full.data$id <- full.data$Animal.ID
full.data$movi_qpcr <- full.data$WADDL.Movi.PCR
full.data$cap_bioyr <- full.data$CAP_BIOYR
full.data <- subset(full.data, movi_qpcr != "Indeterminate")
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")

full.data <- subset(full.data, Age.Class != "Lamb" & Age.Class != "Yearling")

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
data <- subset(full.data , !(id %in% observed.once | id %in% observed.twice) & MostLikelyAge >= 2)
data$id <- factor(data$id)

pcr.by.id.full <- table(full.data$qPCRResult, full.data$id)
prop.pos.full <- prop.table(pcr.by.id.full, margin = 2)[1, ]

pcr.by.id <- table(data$qPCRResult, data$id)
prop.pos <- prop.table(pcr.by.id, margin = 2)[1, ]
rbind(pcr.by.id, 1-prop.pos)


# subset down to adult ewes with all positive results over three or more tests. 
all.pos <- names(prop.pos)[which(prop.pos == 1)]
  # 6 animals: 04Lo49, 05LO53,06LO61, 11LO40, 11LO44, 13LO02

# what were the durations of their infections?
infection.times <- rep(NA, 6)
for(i in 1:6){
  k <- subset(full.data, Animal.ID == as.character(all.pos)[i])
  infection.times[i] <- difftime(as.Date(as.character(k$Capture.Date[dim(k)[1]]), format = "%d/%m/%Y"),
                                 as.Date(as.character(k$Capture.Date[1]), format = "%d/%m/%Y"))
}
mean(infection.times)

# bring in Movi survival data
movi <- read.csv("./Data/MoviSurvivalDate_20160220.csv", header = T)
movi$stopdate.Date <- as.Date(movi$stopdate, format = "%m/%d/%Y")
movi$startdate.Date <- as.Date(movi$startdate, format = "%m/%d/%Y")
movi$previousstartsamp.Date <- as.Date(movi$previousstartsamp, format = "%m/%d/%Y")
movi$previousendsamp.Date <- as.Date(movi$previousendsamp, format = "%m/%d/%Y")

movi$start.mid <- movi$previousstartsamp.Date + floor(movi$startdate.Date - movi$previousstartsamp.Date)/2
movi$start.mid[is.na(movi$previousstartsamp.Date) == T] <- movi$startdate.Date[is.na(movi$previousstartsamp.Date) == T]
movi$stop.mid <- movi$previousendsamp.Date + floor(movi$stopdate.Date - movi$previousendsamp.Date)/2
movi$stop.mid[is.na(movi$previousendsamp.Date) == T] <- movi$stopdate.Date[is.na(movi$previousendsamp.Date) == T]
movi$diffmid <- difftime(movi$stop.mid, movi$start.mid)

which(movi$diffmid >= mean(infection.times))



# build Movi survival data
individs <- vector("list", length(levels(factor(full.data$id))))
infection.statuses <- get.infected <- recover <- last.negative.date <- first.positive.date <- last.positive.date <- first.negative.date <- vector("list", length(individs))
for(i in 1:length(individs)){
  k <- subset(full.data, as.character(id) == as.character(levels(factor(full.data$id)))[i])
  infection.statuses[[i]] <- paste(k$qPCRResult, sep = "", collapse = "")
  get.infected[[i]] <- gregexpr("NP", infection.statuses[[i]])
  last.negative.date[[i]] <- first.positive.date[[i]] <- rep(NA, length(get.infected[[i]][[1]]))
  if(get.infected[[i]][[1]][1] >= 1){
    for(j in 1:length(last.negative.date[[i]])){
      last.negative.date[[i]][j] <- as.character(k$Capture.Date[get.infected[[i]][[1]][j]])
      first.positive.date[[i]][j] <- as.character(k$Capture.Date[get.infected[[i]][[1]][j] + 1])
    }
  }
  
  recover[[i]] <- gregexpr("PN", infection.statuses[[i]])
  last.positive.date[[i]] <- first.negative.date[[i]] <- rep(NA, length(recover[[i]][[1]]))
  if(recover[[i]][[1]][1] >= 1){
    for(j in 1:length(last.positive.date[[i]])){
      last.positive.date[[i]][j] <- as.character(k$Capture.Date[recover[[i]][[1]][j]])
      first.negative.date[[i]][j] <- as.character(k$Capture.Date[recover[[i]][[1]][j] + 1])
    }
  }
}
