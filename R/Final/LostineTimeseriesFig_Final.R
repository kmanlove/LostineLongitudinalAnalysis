compd.data <- read.csv("./Data/compiled_data_summary_151115b.csv", header = T)
los <- subset(compd.data, Pop == "Lostine")
los$CLASS[43:44] <- "LAMBS"
# Lostine population estimate timeseries
pt.col <- ifelse(is.na(los$CLASS) == T, "white", 
                 ifelse(los$CLASS == "HEALTHY", "black", "red"))

# svg("./Plots/LostineTimeSeries_20160622.svg", width = 6, height = 5)
plot(los$PopEst ~ los$year, type = "b", ylab = "Population size", 
     xlab = "Year", las = 1, ylim = c(0, 125))
points(los$PopEst ~ los$year, col = pt.col, pch = 16)
points(los$PopEst ~ los$year, col = "black", pch = 1)
shp.rem <- subset (los, NoSheepRem != 0 & is.na(NoSheepRem) != T)
shp.rel <- subset (los, NoSheepRel != 0 & is.na(NoSheepRel) != T)
text(x = shp.rem$year, y = shp.rem$PopEst + 8, shp.rem$NoSheepRem,
     cex = .60)
#points(los$Lambs ~ los$year, pch = 5)
points(los$Lambs ~ los$year, pch = 17, col = "grey60")
lines(los$Lambs ~ los$year)
leg.text <- c("Infected or suspected infected", "Suspected healthy", "# Recruits", "No health status data")
legend("topright", leg.text, col = c("red", "black", "grey60", "black"), 
       pch = c(16, 16, 17, 1), pt.cex = c(1, 1, 1, 1), bty = "n", cex = .65)
# dev.off()


require(shape)
#svg("./Plots/LostineTimeSeries_20160628_V2.svg", width = 6, height = 5)
plot(los$PopEst ~ los$year, type = "l", ylab = "# animals", 
     xlab = "Year", las = 1, ylim = c(0, 135))
points(los$PopEst ~ los$year, col = pt.col, pch = 16)
points(los$PopEst ~ los$year, col = "black", pch = 1)
shp.rem <- subset (los, NoSheepRem != 0 & is.na(NoSheepRem) != T)
shp.rel <- subset (los, NoSheepRel != 0 & is.na(NoSheepRel) != T)
#text(x = shp.rem$year, y = shp.rem$PopEst + 8, shp.rem$NoSheepRem,
#     cex = .60)
text(x = shp.rem$year[-1], y = shp.rem$PopEst[-1] + 8, shp.rem$NoSheepRem[-dim(shp.rem)[1]],
     cex = .60)
points(shp.rem$PopEst[-1] + shp.rem$NoSheepRem[-dim(shp.rem)[1]] ~ 
         (shp.rem$year[-1]), col = rgb(.5, .5, .5, .5))
points(shp.rem$PopEst[-1] + shp.rem$NoSheepRem[-dim(shp.rem)[1]] ~ 
         (shp.rem$year[-1]), col = rgb(.5, .5, .5, .5))

#points(los$Lambs ~ los$year, pch = 5)
points(los$Lambs ~ los$year, pch = 17, col = "grey60")
lines(los$Lambs ~ los$year)
leg.text <- c("Infected or suspected infected", "Suspected healthy", "Lambs recruited", "No health status data")
legend("topright", leg.text, col = c("red", "black", "grey60", "black"), 
       pch = c(16, 16, 17, 1), pt.cex = c(1, 1, 1, 1), bty = "n", cex = .65)
# require(graphics)
for(i in 2:dim(shp.rem)[1]){
  Arrows(x0 = shp.rem$year[i], x1 = shp.rem$year[i], 
           y0 = (shp.rem$PopEst[i] + shp.rem$NoSheepRem[i-1]),
           y1 = (shp.rem$PopEst[i])+4, cex = .5, code = 2,
         col = rgb(.5, .5, .5, .3))
}
#dev.off()




# # Prop of pos results within-host
# # read in Lostine data
# data <- read.csv("./Data/LostineMoviPastLung_150630_rp.csv", header = T)
# 
# #-- build and explore determination cohort --#
# Animal.tab <- table(data$id)
# #AnimalByOutcome.tab <- table(data$id, data$WADDL.Movi.PCR)
# AnimalByOutcome.tab <- table(data$id, data$movi_qpcr)
# Animal.tab[which(Animal.tab >= 4)]
# Animal.tab[which(Animal.tab == 3)]
# 
# ManySamplesOutcomes.tab <- AnimalByOutcome.tab[which(Animal.tab >= 4), ]
# ManySamplesProps <- ManySamplesOutcomes.tab[ ,1] / (ManySamplesOutcomes.tab[ ,1] + ManySamplesOutcomes.tab[ ,3])
# 
# # histogram of prop positive
# data$new.date <- strptime(data$capture_date, forma = "%m/%d/%Y")
# plot(x = 0, y = 0, cex = 0, ylim = c(0, length(levels(factor(data$id)))))
# for(i in 1:length(levels(factor(data$id)))){
#   k <- subset(data, id == levels(factor(data$id))[i])
#   points(y = rep(i, dim(k)[1]), x = k$new.date, pch = 16, col = as.numeric(k$movi_qpcr))
# }
# 
# cols.in <- c("red", "pink", "grey50")
# plot(data$new.date, as.numeric(data$id), type = "p", 
#      xaxt = "n", pch = 16, 
#      col = cols.in[as.numeric(as.factor(data$movi_qpcr))], 
#      xlab = "date", ylab = "Animal")
# r <- as.POSIXct(round(range(data$new.date), "days"))
# axis.Date(1, at = seq(r[1], r[2], by = "day"), format = "%H")
# 
# par(mfrow = c(1, 1), las = 1, mar = c(4, 4, 2, 2), oma = c(1, 1, 1, 1))
# hist(ManySamplesProps, breaks = 10, col = "grey80", 
#      xlab = "Proportion of PCR-positive tests for \n the 22 animals tested 4 or more times", main = "")
# abline(v = 0.67, col = "red", lty = 2, lwd = 2)
# text(x = 0.3, y = 4.5, "Consistently \n negative")
# text(x = 0.8, y = 4.5, "Consistently \n positive")
