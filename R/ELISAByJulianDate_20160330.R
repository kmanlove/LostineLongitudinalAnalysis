data <- read.csv("./Data/LostineMoviPastLung_150630.csv", header = T)

data$Movipneumonia.ELISA <- as.numeric(as.character(data$Movipneumonia.ELISA))
data$JulianDate <- as.numeric(as.character(data$JulianDate))

plot(data$Movipneumonia.ELISA ~ data$JulianDate)
lines(lowess(data$Movipneumonia.ELISA ~ data$JulianDate, 
             f = 5/7, delta = 0.2 * diff(range(data$JulianDate))), 
      lty = 1, lwd = 2)


