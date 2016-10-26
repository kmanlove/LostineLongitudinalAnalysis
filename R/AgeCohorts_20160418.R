# age-structure likelihood

# basic logic:
# each animal's age is treated as a draw from a multinomial likelihood
# combine info on recr in a given cohort with my age-specific survival estimates
# constrain state space to just be _possible_ birth cohorts for each animal. 
# constrain most likely cohort based on cohorts of known-age animals, 
#   and likelihoods that other animals get assigned to that cohort

compd.data <- read.csv("./Data/compiled_data_summary_151115b.csv", header = T)
los <- subset(compd.data, Pop == "Lostine" & year >= 1990)
# use the posterior survival probs from the BHS demography paper
krm.surv.probs.he <- c(rep(.99, 2),
                       rep(.90, 4),
                       rep(.84, 6),
                       rep(.85, 7),
                       rep(0, 8) 
                       # fill in 0's to be sure vector is long-enough for  next for loop. 
                       ) 

cohort <- matrix(NA, nrow = 27, ncol = 27)
for(i in 1:27){
  k <- subset(los, year == (1989 + i))
    for(j in 1:27){
      cohort[i, j] <- ifelse(j<i, 0, k$Lambs * prod(krm.surv.probs.he[1:(j - i)]))
  }
}



# likelihood an animal comes from a given cohort
samples <- read.csv("./Data/LostineSampleData_20160505_v1.csv", header = T)
studysheep.in <- read.csv("./Data/Study_sheep_20151215.csv", header = T) 
studysheep.in$KnownAge <- ifelse(!(is.na(studysheep.in$Tooth_Age)) | studysheep.in$AgeatEntry < 4, 1, 0)

animals <- start.cohort <- last.cohort <- omit <- most.likely.cohort <- known.age <- rep(NA, length(levels(samples$Animal.ID)))
cohort.dist <- matrix(0, nrow = length(animals), ncol = 27)
years.list <- vector("list", length(animals))
for(i in 1:length(animals)){
  k <- subset(samples, Animal.ID == levels(samples$Animal.ID)[i])
  studysheep <- subset(studysheep.in, as.character(Animal_ID) == levels(samples$Animal.ID)[i])
  animals[i] <- as.character(levels(samples$Animal.ID)[i])
  known.age[i] <- ifelse(is.na(studysheep$Tooth_Age[1]) == T & k$AgeatEntry[1] >= 4, 0, 1)
  if(k$AgeatEntry[1] <= 3.9){
    start.cohort[i] <- floor(k$ENTRY_BIOYR[1] - k$AgeatEntry[1])
    last.cohort[i] <- ceiling(k$ENTRY_BIOYR[1] - k$AgeatEntry[1])
  } else{
    start.cohort[i] <- floor(max(k$CAP_BIOYR) - 20)
    last.cohort[i] <- ceiling(k$ENTRY_BIOYR[1] - k$AgeatEntry[1])
  }
  years.list[[i]] <- seq(floor(start.cohort[i]), ceiling(last.cohort[i]), by = 1)
  omit[i] <- ifelse(max(k$CAP_BIOYR) <= 2005, 1, 0)
  
  if(omit[i] == 0){
    cohort.dist[i, years.list[[i]] - 1989] <- cohort[(years.list[[i]] - 1989), (k$CAP_BIOYR[1] - 1989)] / sum(cohort[(years.list[[i]] - 1989), (k$CAP_BIOYR[1] - 1989)])
  }
  most.likely.cohort[i] <- seq(1990, 2016, by = 1)[which.max(cohort.dist[i, ])]
  print(i)
}

# build animal dataframe with ages
animal.data <- as.data.frame(cbind(as.character(animals), start.cohort, last.cohort, most.likely.cohort, omit, known.age, round(cohort.dist, 3)))
names(animal.data) <- c("Animal.ID", "FirstPossibleCohort", "LastPossibleCohort", "MostLikelyCohort", "Omit", "KnownAge", paste("P.", as.character(seq(1990, 2016, by = 1)), sep= ""))
head(animal.data)

animal.data <- subset(animal.data, Omit == 0)
animal.data$MostLikelyCohortOrig <- animal.data$MostLikelyCohort

# now, adjust for possible number of animals in each cohort
most.likely.cohort.adj <- rep(NA, length(levels(factor(animal.data$Animal.ID))))
for(i in 1:length(levels(factor(animal.data$Animal.ID)))){
  # 1. subset animal data to get all animals assigned to this animal's most likely cohort.
  k <- subset(animal.data, MostLikelyCohort == animal.data$MostLikelyCohort[i])
  animal.likelihoods <- subset(animal.data, Animal.ID == levels(factor(animal.data$Animal.ID))[i])[, 7:33]
  # 2. determine how many of those animals have known ages
  known.age.in.cohort <- dim(subset(k, KnownAge == 1))[1]
  # 3. get the max possible cohort size, using the Lambs field in the lostine compiled data
  possible.animals <- subset(los, year == animal.data$MostLikelyCohort[i])$Lambs
  # 4. figure out how many possible ewes of unknown age exist
  possible.ewes.unknown.age <- ceiling((possible.animals - known.age.in.cohort) / 2)
  # 5. if there are more possible animals of unknown age than there are animals of unknown ages whose 
  #    most likely cohort is this one, assign all those animals to this cohort.
  if(possible.ewes.unknown.age >= (dim(k)[1] - known.age.in.cohort)){
    most.likely.cohort.adj[i] <- as.numeric(as.character(animal.data$MostLikelyCohort[i]))
  } else {
    # 6. pull out likelihoods for each animal in k in the cohort year
    cohort.name <- paste("P.", k$MostLikelyCohort[1], sep = "")
    cohort.yr.likelihoods <- subset(k, select = cohort.name)
    likelihood.ranks <- rank((1 - as.numeric(as.character(cohort.yr.likelihoods[ ,1]))))
    in.this.cohort <- k[likelihood.ranks <= possible.ewes.unknown.age, ]
    if(levels(factor(animal.data$Animal.ID))[i] %in% in.this.cohort$Animal.ID){
      most.likely.cohort.adj[i] <- as.numeric(as.character(k$MostLikelyCohort[1]))
    } else {
      next.highest <- strsplit(x = names(animal.likelihoods)[which.max(as.numeric(as.character(as.vector(t(animal.likelihoods))))[-which.max(as.numeric(as.character(as.vector(t(animal.likelihoods)))))])], split = "P.")[[1]][2]
      animal.data$MostLikelyCohort[which(animal.data$Animal.ID == levels(factor(animal.data$Animal.ID))[i])] <- as.numeric(as.character(next.highest))
    }                                        
#    next.most.likely[i]
    print(i)
  }
}
table(is.na(most.likely.cohort.adj))
# iterate through the above for-loop until there are no NAs left. 

# write out the csv file (this is the csv I sent around)
# write.csv(animal.data, "./Data/AnimalAgeData_20160521_krm.csv")


# plot out the age distributions to be sure MostLikelyCohortOrig and Most LikelyCohort aren't identical
par(mfrow = c(1, 2))
hist(as.numeric(as.character(animal.data$MostLikelyCohortOrig)), breaks = 20)
hist(as.numeric(as.character(animal.data$MostLikelyCohort)), breaks = 20)

# attach animal cohort data to animal handling dataframe (to associate ages with each sampling event)
full.data$MostLikelyAge <- rep(NA, dim(full.data)[1])
for(i in 1:dim(full.data)[1]){
  k <- subset(animal.data, as.character(Animal.ID) == as.character(full.data$Animal.ID)[i])
  full.data$MostLikelyAge[i] <- as.numeric(as.character(full.data$CAP_BIOYR[i])) - as.numeric(as.character(k$MostLikelyCohort[1]))
  print(i)
}

#write.csv(full.data, "./Data/LostineSampleData_20160505_v1.csv")


# KRM age structure plots and analysis of quadrature (DON'T USE MY ANALYSIS OF QUADRATURE; DAVID'S/ALAN'S IS BETTER)
rainbow.cols <- rainbow(n = 27)

par(mfrow = c(1, 3), las = 1, cex.axis = 1.2, cex.lab = 1.5, mar = c(5, 5, 2, 2))
plot(cohort[1, ] ~ seq(1990:2016), type = "l", ylim = c(0, 25), 
     col = rainbow.cols[1], 
     lwd = 2, xaxt = "n",
     main = "Since 1990",
     xlab = "Sampling year", 
     ylab = "Expected number of animals from cohort")
axis(side = 1, at = seq(1990:2016), labels = c(1990:2016))
for(i in 2:27){
  lines(cohort[i, i:27] ~ seq(1990:2016)[-c(1:(i - 1))], 
        col = rainbow.cols[i],
        lwd = 2)
  x.label <- (1989 + i)
  text(x = i, y = cohort[i, i] + .5, labels = as.character(x.label), 
       col = rainbow.cols[i])
}

plot(cohort[1, 18:27] ~ seq(2007,2016, by = 1), type = "l", ylim = c(0, 25), 
     col = rainbow.cols[1], 
     lwd = 2, xaxt = "n",
     main = "Last 10 years",
     xlab = "Sampling year", 
     ylab = "Expected number of animals from cohort")
axis(side = 1, at = seq(2007,2016, by = 1), labels = c(2007:2016))
for(i in 2:27){
  start <- ifelse(i >= 18, i, 18)
  #  seq.in <- ifelse(i >= 18, seq(2007:2016)[-c(1:(i - 1))], seq())
  lines(cohort[i, start:27] ~ seq(1990,2016, by = 1)[-c(1:(start - 1))], 
        col = rainbow.cols[i],
        lwd = 2)
  x.label <- (1989 + i)
  text(x = (1989 + start), y = cohort[i, start] + .5, labels = as.character(x.label), 
       col = rainbow.cols[i])
}

known.age.sheep <- subset(animal.data, KnownAge == 1)
hist.breaks <- hist(as.numeric(as.character(animal.data$MostLikelyCohort)), 
                    col = rgb(0, 0, 0, alpha = .5), #col = rainbow.cols, 
                    breaks = 20, xlab = "Most likely birth year",
                    main = "Lostine sampled animals")$breaks

hist(as.numeric(as.character(known.age.sheep$MostLikelyCohort)), breaks = 20)
par(mfrow = c(1, 1), oma = c(2, 2, 0, 0))
hist(as.numeric(as.character(animal.data$MostLikelyCohort)), 
     col = rgb(.1, .1, .1, alpha = .35), #col = rainbow.cols, 
     breaks = hist.breaks, xlab = "Most likely birth year",
      ylim = c(0, 20), main = "")
par(new = T)
hist(as.numeric(as.character(known.age.sheep$MostLikelyCohort)), 
     col = "grey30", #col = rainbow.cols, 
     breaks = hist.breaks, xlab = "Most likely birth year",
     ylim = c(0, 20), main = "")
leg.text <- c("Known age", "Estimated age")
legend("topleft", leg.text, fill = c("grey30", rgb(.1, .1, .1, alpha = .35)), bty = "n")

hist(as.numeric(as.character(full.data$MostLikelyAge)), 
     col = "grey60", #col = rainbow.cols, 
     breaks = 20, xlab = "Most likely age",
     main = "Lostine sampled animals", xlim = c(0, 20))

# append most likely age to full.data
samples$animal.current.age <- rep(NA, dim(samples)[1])
for(i in 1:dim(samples)[1]){
  k <- subset(animal.data, as.character(Animal.ID) == as.character(samples$Animal.ID)[i])
  samples$animal.current.age[i] <- as.numeric(as.character(samples$CAP_BIOYR[i])) - as.numeric(as.character(k$MostLikelyCohort))
  print(i)
}

samples$movi.pcr.pos <- ifelse(samples$WADDL.Movi.PCR == "Detected" | samples$WADDL.Movi.PCR == "Indeterminate", 1,
                        ifelse(samples$WADDL.Movi.PCR == "Not detected", 0, NA))

samples$animal.current.age.2yrbins <- ifelse(samples$animal.current.age %in% c(0, 1), .5,
                                     ifelse(samples$animal.current.age %in% c(2, 3), 2.5,
                                     ifelse(samples$animal.current.age %in% c(4, 5), 4.5,
                                     ifelse(samples$animal.current.age %in% c(6, 7), 6.5, 
                                     ifelse(samples$animal.current.age %in% c(8, 9), 8.5,
                                     ifelse(samples$animal.current.age %in% c(10, 11), 10.5, 
                                     ifelse(samples$animal.current.age %in% c(12, 13), 12.5,
                                     ifelse(samples$animal.current.age %in% c(14, 15), 14.5, 
                                     ifelse(samples$animal.current.age %in% c(16, 17), 16.5,
                                     ifelse(samples$animal.current.age %in% c(18, 19), 18.5,
                                     ifelse(samples$animal.current.age %in% c(20), 20, NA)))))))))))
require(lme4)
require(graphics)
logist.cat <- glm(movi.pcr.pos ~ factor(animal.current.age), 
                     family = "binomial", data = samples)
logist.cat.2yr <- glm(movi.pcr.pos ~ factor(animal.current.age.2yrbins), 
                  family = "binomial", data = samples)
newdata <- data.frame(animal.current.age.2yrbins = as.numeric(as.character(levels(factor(samples$animal.current.age.2yrbins)))))
logist.cat.preds <- predict(logist.cat.2yr, newdata, type = "response", se = T)
logist.cat.preds <- predict(logist.cat.2yr, newdata, type = "response", se = T)
logist.quad <- glm(movi.pcr.pos ~ animal.current.age + I(animal.current.age^2), 
                   family = "binomial", data = samples)
newdata.quad <- data.frame(animal.current.age = seq(0, 20, length.out = 200))
logist.quad.pred <- predict(logist.quad, newdata.quad, type = "response", se = T)
prev.by.age <- table(samples$movi.pcr.pos, samples$animal.current.age.2yrbins)[2, ]/table(samples$animal.current.age.2yrbins)
ages <- as.numeric(as.character(names(table(samples$animal.current.age.2yrbins))))

# svg("./Plots/AgePrevalenceCurve_20160505.svg", height = 5, width = 6)
par(mfrow = c(1,1))
plot(as.numeric(as.character(prev.by.age)) ~ ages, xlim = c(0, 20), ylim = c(0, 1),
     ylab = "M.ovi prevalence", xlab = "Animal age")
for(i in 1:dim(newdata)[1]){
  segments(x0 = newdata$animal.current.age[i], x1 = newdata$animal.current.age[i], 
        y0 = logist.cat.preds$fit[i] - logist.cat.preds$se.fit[i], 
        y1 = logist.cat.preds$fit[i] + logist.cat.preds$se.fit[i])
  segments(x0 = newdata$animal.current.age[i] - .1, x1 = newdata$animal.current.age[i] + .1, 
           y0 = logist.cat.preds$fit[i], y1 = logist.cat.preds$fit[i])
  text(x = newdata$animal.current.age[i], y = 0.0, cex = .8,
       paste("(", table(samples$animal.current.age.2yrbins)[i], ")", sep = ""))
}
lines(logist.quad.pred$fit ~ newdata.quad$animal.current.age, type = "l")
polygon(x = c(newdata.quad$animal.current.age, rev(newdata.quad$animal.current.age)),
        y = c(logist.quad.pred$fit - 1.96 * logist.quad.pred$se.fit, 
              rev(logist.quad.pred$fit + 1.96 * logist.quad.pred$se.fit)), 
        col = rgb(0, 0, 0, alpha = .2), border = rgb(0, 0, 0, alpha = .4))
# dev.off()


