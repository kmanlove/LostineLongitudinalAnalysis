full.data <- read.csv("./Data/LostineSample2016Final_WithClasses_20160929.csv", header = T)
full.data <- subset(full.data, cap_bioyr >= 2011)
full.data <- subset(full.data, movi_qpcr %in% c("Detected", "Not detected"))
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")

global.min.date <- min(strptime(as.character(full.data$capture_date), format = "%m/%d/%Y"))
full.data$days_into_study <- as.numeric(difftime(strptime(as.character(full.data$capture_date), format = "%m/%d/%Y"),
                                      global.min.date,
                                      units = "days"
                                      ))
         
not.adults <- subset(full.data, age_class != "Adult")
adults <- subset(full.data, age_class == "Adult" & !(id %in% levels(factor(not.adults$id))))
n.adults <- length(levels(factor(adults$id)))
adults$id <- factor(adults$id)

# get animal's classification as neg/inter/carrier on two consec years definition
animal.class <- rep(NA, length(levels(factor(adults$id))))
for(i in 1:length(animal.class)){
  k <- subset(adults, id == levels(factor(adults$id))[i])
  animal.class[i] <- as.character(k$carrier_posConsecYears[1])
}

animal.data <- as.data.frame(cbind(levels(factor(adults$id)), 
                                   as.character(animal.class)))
names(animal.data) <- c("id", "class")

# pull out ids in each group (carriers, inters, negs)
carriers <- subset(animal.data, class == "carrier")
inters <- subset(animal.data, class == "intermittent")
negs <- subset(animal.data, class == "negative")
NAs <- subset(animal.data, is.na(class) == T)

# re-order animal IDs so that they're blocked by carriage status. 
adults$id.2 <- factor(adults$id, 
                      levels = c(levels(factor(NAs$id)),
                                 levels(factor(negs$id)),
                                 levels(factor(inters$id)),
                                 levels(factor(carriers$id))
                      ))

# get time deviations since first positive test for each animal. 
adults$time.rel.to.1st.pos <- rep(NA, dim(adults)[1])
for(i in 1:dim(adults)[1]){
  k <- subset(adults, id == adults$id[i] & FirstPosInd == 1)
  if(dim(k)[1] == 1){
    adults$time.rel.to.1st.pos[i] <- difftime(time2 = strptime(k$capture_date[1], format = "%m/%d/%Y"),
                                              time1 = strptime(adults$capture_date[i], format = "%m/%d/%Y"))
    print("k")
    } else {
    m <- subset(adults, id == adults$id[i] & FirstTestInd == 1)
    adults$time.rel.to.1st.pos[i] <- difftime(time2 = strptime(m$capture_date[1], format = "%m/%d/%Y"),
                                              time1 = strptime(adults$capture_date[i], format = "%m/%d/%Y"))
    
    print("m")
  }
  print(adults$time.rel.to.1st.pos[i])
}

#-------------------------------------#
#-- Sampling intensities for adults --#
#-------------------------------------#

#svg("./Plots/SamplingIntensityFig_20161006_V2.svg", height = 7, width = 6)
par(mfrow = c(1, 2), las = 2, mar = c(6, 6, 1, 1), cex.axis = .7)

# first panel: samples within individuals through time. 
plot(x = 0, y = 0, xlim = c(0, 1700), ylim = c(1, n.adults), cex = 0,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:n.adults){
    abline(h = i, lty = 1, lwd = 1, col = "grey60")
}
points(adults$id.2 ~ adults$days_into_study, pch = 16, 
       col = ifelse(adults$movi_qpcr == "Detected", "red", "black"))
axis(side = 2, at = seq(1, n.adults), labels = levels(factor(adults$id.2)), cex = .4)
axis(side = 1, at = c(260, 260 + 365, 260 + 365*2, 260 + 365*3),
     labels = c("Jan 1, 2013", "Jan 1, 2014", "Jan 1, 2015", "Jan 1, 2016"))

# second panel: first pos = 0
plot(x = 0, y = 0, xlim = c(-365*3, 365*3.2), ylim = c(1, n.adults), cex = 0,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:n.adults){
  abline(h = i, lty = 1, lwd = 1, col = "grey60")
}
points(adults$id.2 ~ adults$time.rel.to.1st.pos, pch = 16, 
       col = ifelse(adults$movi_qpcr == "Detected", "red", "black"))
axis(side = 2, at = seq(1, n.adults), labels = levels(factor(adults$id.2)), cex = .4)
axis(side = 1, at = c(-365*3, -365*2, -365, 0, 365, 365*2, 365*3),
     labels = c(-3, -2, -1, 0, 1, 2, 3))
#dev.off()


#----------------------------------------#
#-- Sampling intensities for juveniles --#
#----------------------------------------#

# get animal's classification as neg/inter/carrier on two consec years definition
juvs <- subset(full.data, !(id %in% levels(factor(adults$id))))
n.juvs <- length(levels(factor(juvs$id)))
juv.class <- rep(NA, length(levels(factor(juvs$id))))
for(i in 1:length(juv.class)){
  k <- subset(juvs, id == levels(factor(juvs$id))[i])
  juv.class[i] <- as.character(k$carrier_posConsecYears[1])
}

juv.data <- as.data.frame(cbind(levels(factor(juvs$id)), 
                                   as.character(juv.class)))
names(juv.data) <- c("id", "class")

# pull out ids in each group (carriers, inters, negs)
juv.carriers <- subset(juv.data, class == "carrier")
juv.inters <- subset(juv.data, class == "intermittent")
juv.negs <- subset(juv.data, class == "negative")
juv.NAs <- subset(juv.data, is.na(class) == T)

# re-order animal IDs so that they're blocked by carriage status. 
juvs$id.2 <- factor(juvs$id, 
                      levels = c(levels(factor(juv.NAs$id)),
                                 levels(factor(juv.negs$id)),
                                 levels(factor(juv.inters$id)),
                                 levels(factor(juv.carriers$id))
                      ))

# get time deviations since first positive test for each animal. 
juvs$time.rel.to.1st.pos <- rep(NA, dim(juvs)[1])
for(i in 1:dim(juvs)[1]){
  k <- subset(juvs, id == juvs$id[i] & FirstPosInd == 1)
  if(dim(k)[1] == 1){
    juvs$time.rel.to.1st.pos[i] <- difftime(time2 = strptime(k$capture_date[1], format = "%m/%d/%Y"),
                                              time1 = strptime(juvs$capture_date[i], format = "%m/%d/%Y"))
    print("k")
  } else {
    m <- subset(juvs, id == juvs$id[i] & FirstTestInd == 1)
    juvs$time.rel.to.1st.pos[i] <- difftime(time2 = strptime(m$capture_date[1], format = "%m/%d/%Y"),
                                              time1 = strptime(juvs$capture_date[i], format = "%m/%d/%Y"))
    
    print("m")
  }
  print(juvs$time.rel.to.1st.pos[i])
}

# get time since birth for each animal. 
juvs$age_days <- rep(NA, dim(juvs)[1])
for(i in 1:dim(juvs)[1]){
    juvs$age_days[i] <- difftime(time2 = strptime(juvs$date_birth[i], format = "%m/%d/%Y"),
                                            time1 = strptime(juvs$capture_date[i], format = "%m/%d/%Y"))
}  
  
#svg("./Plots/SamplingIntensityFig_JUVS_20161006_V2.svg", height = 5.5, width = 6)
par(mfrow = c(1, 2), las = 2, mar = c(6, 6, 1, 1), cex.axis = .7)

# first panel: samples within individuals through time. 
plot(x = 0, y = 0, xlim = c(0, 1700), ylim = c(1, n.juvs), cex = 0,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:n.juvs){
  abline(h = i, lty = 1, lwd = 1, col = "grey60")
}
points(juvs$id.2 ~ juvs$days_into_study, pch = 16, 
       col = ifelse(juvs$movi_qpcr == "Detected", "red", "black"))
axis(side = 2, at = seq(1, n.juvs), labels = levels(factor(juvs$id.2)), cex = .4)
axis(side = 1, at = c(260, 260 + 365, 260 + 365*2, 260 + 365*3),
     labels = c("Jan 1, 2013", "Jan 1, 2014", "Jan 1, 2015", "Jan 1, 2016"))

# second panel: first pos = 0
plot(x = 0, y = 0, xlim = c(-365*3, 365*3.2), ylim = c(1, n.juvs), cex = 0,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:n.juvs){
  abline(h = i, lty = 1, lwd = 1, col = "grey60")
}
points(juvs$id.2 ~ juvs$time.rel.to.1st.pos, pch = 16, 
       col = ifelse(juvs$movi_qpcr == "Detected", "red", "black"))
axis(side = 2, at = seq(1, n.juvs), labels = levels(factor(juvs$id.2)), cex = .4)
axis(side = 1, at = c(-365*3, -365*2, -365, 0, 365, 365*2, 365*3),
     labels = c(-3, -2, -1, 0, 1, 2, 3))
#dev.off()


#---------------------------------#
#-- Juveniles with three panels --#
#---------------------------------#

#svg("./Plots/SamplingIntensityFig_3PanelJUVS_20161006_V1.svg", height = 6, width = 9)
par(mfrow = c(1, 3), las = 2, mar = c(6, 6, 1, 1), cex.axis = .7)

# first panel: samples within individuals through time. 
plot(x = 0, y = 0, xlim = c(0, 1700), ylim = c(1, n.juvs), cex = 0,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:n.juvs){
  abline(h = i, lty = 1, lwd = 1, col = "grey60")
}
points(juvs$id.2 ~ juvs$days_into_study, pch = 16, 
       col = ifelse(juvs$movi_qpcr == "Detected", "red", "black"), cex = 1.5)
axis(side = 2, at = seq(1, n.juvs), labels = levels(factor(juvs$id.2)))
axis(side = 1, at = c(260, 260 + 365, 260 + 365*2, 260 + 365*3),
     labels = c("Jan 1, 2013", "Jan 1, 2014", "Jan 1, 2015", "Jan 1, 2016"))

# second panel: first pos = 0
plot(x = 0, y = 0, xlim = c(-365*3, 365*3.2), ylim = c(1, n.juvs), cex = 0,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:n.juvs){
  abline(h = i, lty = 1, lwd = 1, col = "grey60")
}
points(juvs$id.2 ~ juvs$time.rel.to.1st.pos, pch = 16, 
       col = ifelse(juvs$movi_qpcr == "Detected", "red", "black"), cex = 1.5)
axis(side = 2, at = seq(1, n.juvs), labels = levels(factor(juvs$id.2)))
axis(side = 1, at = c(-365*3, -365*2, -365, 0, 365, 365*2, 365*3),
     labels = c(-3, -2, -1, 0, 1, 2, 3))

# third panel: birth day = 0
plot(x = 0, y = 0, xlim = c(0, 365*4), ylim = c(1, n.juvs), cex = 0,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:n.juvs){
  abline(h = i, lty = 1, lwd = 1, col = "grey60")
}
points(juvs$id.2 ~ juvs$age_days, pch = 16, 
       col = ifelse(juvs$movi_qpcr == "Detected", "red", "black"), cex = 1.5)
axis(side = 2, at = seq(1, n.juvs), labels = levels(factor(juvs$id.2)))
axis(side = 1, at = c(0, 365, 365*2, 365*3, 365*4),
     labels = c( 0, 1, 2, 3, 4))
#dev.off()