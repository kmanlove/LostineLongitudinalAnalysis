full.data <- read.csv("./Data/LostineSample2016Final_WithClasses_20160927.csv", header = T)
ind.data <- subset(full.data, movi_qpcr == "Indeterminate")
full.data <- subset(full.data, movi_qpcr %in% c("Detected", "Not detected"))
full.data$qPCRResult <- ifelse(full.data$movi_qpcr == "Detected" | full.data$movi_qpcr == "Indeterminate", "P", "N")
full.data <- subset(full.data, cap_bioyr >= 2008)

# remove all animals with only one observation
observed.once <- names(table(full.data$id))[which(table(full.data$id) == 1)]

data <- subset(full.data, !(id %in% observed.once))
data$id <- factor(data$id)
data$movi_qpcr <- ifelse(data$movi_qpcr == "Indeterminate", "Detected", as.character(data$movi_qpcr))
# get string lengths of all consecutive same outcomes
individ.list <- vector("list", length(levels(factor(data$id))))
string.length <- string.status <- vector("list", length(levels(data$id)))
for(i in 1:length(levels(data$id))){
  individ.list[[i]] <- subset(data, id == levels(data$id)[i])
  switches <- which(individ.list[[i]]$movi_qpcr[-dim(individ.list[[i]])[1]] != individ.list[[i]]$movi_qpcr[-1]) + 1
  switches.full <- c(1, switches, dim(individ.list[[i]])[1])
  string.length[[i]] <- string.status[[i]] <- rep(NA, length(switches.full) - 1)
  if(length(switches.full) == 1){
    string.length[[i]][1] <- difftime(time1 = strptime(individ.list[[i]]$capture_date[dim(individ.list[[i]])[1]], format = "%m/%d/%Y"),
                                      time2 = strptime(individ.list[[i]]$capture_date[1], format = "%m/%d/%Y")
    )
    string.status[[i]][1] <- individ.list[[i]]$movi_qpcr[1]
    
  } else {
    for(j in 1:(length(switches.full) - 1)){
          string.length[[i]][j] <- difftime(time1 = strptime(individ.list[[i]]$capture_date[switches.full[j + 1]], format = "%m/%d/%Y"),
                                            time2 = strptime(individ.list[[i]]$capture_date[switches.full[j]], format = "%m/%d/%Y")
          )
          string.status[[i]][j] <- as.character(individ.list[[i]][switches.full[j], ]$movi_qpcr)
    }
  }
  print(i)
}

full.string.length.orig <- do.call("c", string.length)
full.string.status.orig <- do.call("c", string.status)
orig.mean.length <- mean(full.string.length.orig)
orig.mean.length.detected <- mean(full.string.length.orig[which(full.string.status.orig == "Detected")])
orig.mean.length.undetected <- mean(full.string.length.orig[which(full.string.status.orig == "Not detected")])

# get string lengths of consecutive same outcomes for animals classified as CARRIERS
carrier1 <- subset(data, carrier_posConsecYears == "carrier")
carrier1$id <- factor(carrier1$id)
individ.list.carrier1 <- vector("list", length(levels(factor(carrier1$id))))
string.length.carrier1 <- string.status.carrier1 <- vector("list", length(levels(carrier1$id)))
for(i in 1:length(levels(factor(carrier1$id)))){
  individ.list.carrier1[[i]] <- subset(carrier1, id == levels(carrier1$id)[i])
  switches.carrier1 <- which(individ.list.carrier1[[i]]$movi_qpcr[-dim(individ.list.carrier1[[i]])[1]] != individ.list.carrier1[[i]]$movi_qpcr[-1]) + 1
  switches.full.carrier1 <- c(1, switches.carrier1, dim(individ.list.carrier1[[i]])[1])
  string.length.carrier1[[i]] <- string.status.carrier1[[i]] <- rep(NA, length(switches.full.carrier1) - 1)
  if(length(switches.full.carrier1) == 1){
    string.length.carrier1[[i]][1] <- difftime(time1 = strptime(individ.list.carrier1[[i]]$capture_date[dim(individ.list.carrier1[[i]])[1]], format = "%m/%d/%Y"),
                                      time2 = strptime(individ.list.carrier1[[i]]$capture_date[1], format = "%m/%d/%Y")
    )
    string.status.carrier1[[i]][1] <- individ.list.carrier1[[i]]$movi_qpcr[1]
    
  } else {
    for(j in 1:(length(switches.full.carrier1) - 1)){
      string.length.carrier1[[i]][j] <- difftime(time1 = strptime(individ.list.carrier1[[i]]$capture_date[switches.full.carrier1[j + 1]], format = "%m/%d/%Y"),
                                        time2 = strptime(individ.list.carrier1[[i]]$capture_date[switches.full.carrier1[j]], format = "%m/%d/%Y")
      )
      string.status.carrier1[[i]][j] <- as.character(individ.list.carrier1[[i]][switches.full.carrier1[j], ]$movi_qpcr)
    }
  }
  print(i)
}

full.string.length.orig.carrier1 <- do.call("c", string.length.carrier1)
full.string.status.orig.carrier1 <- do.call("c", string.status.carrier1)
orig.mean.length.carrier1 <- mean(full.string.length.orig.carrier1)
orig.mean.length.detected.carrier1 <- mean(full.string.length.orig.carrier1[which(full.string.status.orig.carrier1 == "Detected")])
orig.mean.length.undetected.carrier1 <- mean(full.string.length.orig.carrier1[which(full.string.status.orig.carrier1 == "Not detected")])


# get string lengths of consecutive same outcomes for animals classified as CARRIERS (pos two consec years)
carrier2 <- subset(data, carrier_posConsecYears == "carrier")
carrier2$id <- factor(carrier2$id)
individ.list.carrier2 <- vector("list", length(levels(factor(carrier2$id))))
string.length.carrier2 <- string.status.carrier2 <- vector("list", length(levels(carrier2$id)))
for(i in 1:length(levels(factor(carrier2$id)))){
  individ.list.carrier2[[i]] <- subset(carrier2, id == levels(carrier2$id)[i])
  switches.carrier2 <- which(individ.list.carrier2[[i]]$movi_qpcr[-dim(individ.list.carrier2[[i]])[1]] != individ.list.carrier2[[i]]$movi_qpcr[-1]) + 1
  switches.full.carrier2 <- c(1, switches.carrier2, dim(individ.list.carrier2[[i]])[1])
  string.length.carrier2[[i]] <- string.status.carrier2[[i]] <- rep(NA, length(switches.full.carrier2) - 1)
  if(length(switches.full.carrier2) == 1){
    string.length.carrier2[[i]][1] <- difftime(time1 = strptime(individ.list.carrier2[[i]]$capture_date[dim(individ.list.carrier2[[i]])[1]], format = "%m/%d/%Y"),
                                               time2 = strptime(individ.list.carrier2[[i]]$capture_date[1], format = "%m/%d/%Y")
    )
    string.status.carrier2[[i]][1] <- individ.list.carrier2[[i]]$movi_qpcr[1]
    
  } else {
    for(j in 1:(length(switches.full.carrier2) - 1)){
      string.length.carrier2[[i]][j] <- difftime(time1 = strptime(individ.list.carrier2[[i]]$capture_date[switches.full.carrier2[j + 1]], format = "%m/%d/%Y"),
                                                 time2 = strptime(individ.list.carrier2[[i]]$capture_date[switches.full.carrier2[j]], format = "%m/%d/%Y")
      )
      string.status.carrier2[[i]][j] <- as.character(individ.list.carrier2[[i]][switches.full.carrier2[j], ]$movi_qpcr)
    }
  }
  print(i)
}

full.string.length.orig.carrier2 <- do.call("c", string.length.carrier2)
full.string.status.orig.carrier2 <- do.call("c", string.status.carrier2)
orig.mean.length.carrier2 <- mean(full.string.length.orig.carrier2)
orig.mean.length.detected.carrier2 <- mean(full.string.length.orig.carrier2[which(full.string.status.orig.carrier2 == "Detected")])
orig.mean.length.undetected.carrier2 <- mean(full.string.length.orig.carrier2[which(full.string.status.orig.carrier2 == "Not detected")])


# generate bootstrapped distribution
data <- subset(data, is.na(movi_qpcr) == F & !(movi_qpcr == ""))
data.boot.fun <- function(nboot, data){
  # build storage objects
  full.string.length <- full.string.status <- data.list <- vector("list", nboot)
  boot.mean.length <- boot.mean.length.detected <- boot.mean.length.undetected <- rep(NA, nboot)
  
  individ.list <- vector("list", length(levels(factor(data$id))))
  string.length <- string.status <- vector("list", length(levels(data$id)))
  
  for(n in 1:nboot){
    data.list[[n]] <- data
    data.list[[n]]$movi_qpcr <- sample(data$movi_qpcr, size = dim(data)[1], replace = F)
  for(i in 1:length(levels(data.list[[n]]$id))){
    individ.list[[i]] <- subset(data.list[[n]], id == levels(data.list[[n]]$id)[i])
    switches <- which(individ.list[[i]]$movi_qpcr[-dim(individ.list[[i]])[1]] != individ.list[[i]]$movi_qpcr[-1]) + 1
    switches.full <- c(1, switches, dim(individ.list[[i]])[1])
    string.length[[i]] <- string.status[[i]] <- rep(NA, length(switches.full) - 1)
    if(length(switches.full) == 1){
      string.length[[i]][1] <- difftime(time1 = strptime(individ.list[[i]]$capture_date[dim(individ.list[[i]])[1]], format = "%m/%d/%Y"),
                                        time2 = strptime(individ.list[[i]]$capture_date[1], format = "%m/%d/%Y")
      )
      string.status[[i]][1] <- individ.list[[i]]$movi_qpcr[1]
      
    } else {
      for(j in 1:(length(switches.full) - 1)){
        string.length[[i]][j] <- difftime(time1 = strptime(individ.list[[i]]$capture_date[switches.full[j + 1]], format = "%m/%d/%Y"),
                                          time2 = strptime(individ.list[[i]]$capture_date[switches.full[j]], format = "%m/%d/%Y")
        )
        string.status[[i]][j] <- as.character(individ.list[[i]][switches.full[j], ]$movi_qpcr)
      }
    }
  }
  full.string.length[[n]] <- do.call("c", string.length)
  full.string.status[[n]] <- do.call("c", string.status)
  
  boot.mean.length[n] <- mean(full.string.length[[n]])
  boot.mean.length.detected[n] <- mean(full.string.length[[n]][which(full.string.status[[n]] == "Detected")])
  boot.mean.length.undetected[n] <- mean(full.string.length[[n]][which(full.string.status[[n]] == "Not detected")])
  
  print(n)
  }
  out.list <- list(full.string.length = full.string.length,
                   full.string.status = full.string.status,
                   boot.mean.length = boot.mean.length,
                   boot.mean.length.detected = boot.mean.length.detected,
                   boot.mean.length.undetected = boot.mean.length.undetected)
  return(out.list)
}

data.boot.test <- data.boot.fun(nboot = 5000, data = data)
table(data.boot.test$boot.mean.length.detected >= orig.mean.length.detected)

#svg("./Plots/StatusPermutationTests_20160624.svg", height = 3, width = 7)
par(mfrow = c(1, 3), las = 1)
hist(data.boot.test$boot.mean.length, 
     breaks = seq(100, 600, length.out = 30), col = "grey80",
     main = "", xlab = "Mean length of same status", ylim = c(0, 2500))
abline(v = orig.mean.length, lwd = 2, col = "red")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n")
hist(data.boot.test$boot.mean.length.detected, 
     breaks = seq(100, 500, length.out = 30), col = "grey80", ylim = c(0, 1000),
     main = "", xlab = "Mean length of time qPCR-positive")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n")
abline(v = orig.mean.length.detected, lwd = 2, col = "red")
hist(data.boot.test$boot.mean.length.undetected, ylim = c(0, 1000),
     breaks = seq(100, 500, length.out = 30), col = "grey80",
     main = "", xlab = "Mean length of time qPCR-negative")
abline(v = orig.mean.length.undetected, lwd = 2, col = "red")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n")
#dev.off()

#svg("./Plots/StatusPermutationTests_20160624.svg", height = 3, width = 7)
par(mfrow = c(2, 2), las = 1)
hist(data.boot.test$boot.mean.length, 
     breaks = seq(100, 600, length.out = 30), col = "grey80",
     main = "", xlab = "Mean length of same status", ylim = c(0, 2500))
abline(v = orig.mean.length, lwd = 2, col = "red")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n")

hist(data.boot.test$boot.mean.length.undetected, ylim = c(0, 1500),
     breaks = seq(100, 600, length.out = 30), col = "grey80",
     main = "", xlab = "Mean length of time qPCR-negative")
abline(v = orig.mean.length.undetected, lwd = 2, col = "red")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n")

hist(data.boot.test$boot.mean.length.detected, 
     breaks = seq(100, 600, length.out = 30), col = "grey80", ylim = c(0, 1500),
     main = "", xlab = "Mean length of time qPCR-positive \n(2 + in one year)")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n", cex = .9)
abline(v = orig.mean.length.detected.carrier1, lwd = 2, col = "red")

hist(data.boot.test$boot.mean.length.detected, 
     breaks = seq(100, 600, length.out = 30), col = "grey80", ylim = c(0, 1500),
     main = "", xlab = "Mean length of time qPCR-positive \n(+ in consec. years)")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n", cex = .9)
abline(v = orig.mean.length.detected.carrier2, lwd = 2, col = "red")

#dev.off()

#svg("./Plots/StatusPermutationTests_20160930.svg", height = 3, width = 7)
par(mfrow = c(1, 3), las = 1)
hist(data.boot.test$boot.mean.length, 
     breaks = seq(100, 600, length.out = 30), col = "grey80",
     main = "", xlab = "Mean length of same status", ylim = c(0, 2500))
abline(v = orig.mean.length, lwd = 2, col = "red")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n", cex = .9)

hist(data.boot.test$boot.mean.length.undetected, ylim = c(0, 1500),
     breaks = seq(100, 600, length.out = 30), col = "grey80",
     main = "", xlab = "Mean length of time qPCR-negative")
abline(v = orig.mean.length.undetected, lwd = 2, col = "red")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n", cex = .9)

hist(data.boot.test$boot.mean.length.detected, 
     breaks = seq(100, 600, length.out = 30), col = "grey80", ylim = c(0, 1500),
     main = "", xlab = "Mean length of time qPCR-positive")
leg.text <- c("Bootstrapped values", "Observed value")
legend("topleft", leg.text, fill = c("grey60", "red"), bty = "n", cex = .9)
#abline(v = orig.mean.length.detected, lwd = 2, col = "red")
abline(v = orig.mean.length.detected.carrier1, lwd = 1, col = "red", lty = 2)
#abline(v = orig.mean.length.detected.carrier2, lwd = 2, col = "red", lty = 3)
leg.text.2 <- c("+ in two consec years")
legend(x = 90, y = 1400, leg.text.2, lty = c(1, 2, 3), lwd = rep(2, 3), col = rep("red", 3), bty = "n", cex = .9)

#dev.off()