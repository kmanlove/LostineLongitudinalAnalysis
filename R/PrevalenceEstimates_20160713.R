# Movi prevalence estimation
samples <- read.csv("./Data/LostineSample2016Final.csv", header = T)
length(levels(factor(samples$id)))

samples.small <- subset(samples, age_class == "Adult" & sex == "F")

# pull off one obs per animal
ind.carriage.class.1 <- ind.carriage.class.2 <- ind.age.class <- ind.id <- rep(NA, length(levels(samples.small$id)))
for(i in 1:length(levels(samples.small$id))){
  k <- subset(samples.small, as.character(id) == as.character(levels(samples.small$id))[i])
  ind.carriage.class.1[i] <- as.character(k$carriage_twoPosOneYear)[1]
  ind.carriage.class.2[i] <- as.character(k$carrier_posConsecYears)[1]
  ind.id[i] <- as.character(levels(samples.small$id))[i]
}

carriage.status.props.1 <- table(na.omit(ind.carriage.class.1))/sum(table(na.omit(ind.carriage.class.1)))
carriage.status.props.2 <- table(na.omit(ind.carriage.class.2))/sum(table(na.omit(ind.carriage.class.2)))

nboot <- 1000
carriage.1.boot <- carriage.2.boot <- matrix(NA, nrow = nboot, ncol = 3)
ind.carriage.class.1.nonas <- na.omit(ind.carriage.class.1)
ind.carriage.class.2.nonas <- na.omit(ind.carriage.class.2)
for(i in 1:nboot){
  k <- sample(1:length(ind.carriage.class.2.nonas), size = length(ind.carriage.class.2.nonas), replace = T)
  carriage.1.boot[i, ] <- table(factor(ind.carriage.class.1.nonas[k], levels = c("carrier", "intermittent", "negative")))/sum(table(factor(ind.carriage.class.1.nonas[k], levels = c("carrier", "intermittent", "negative"))))
  carriage.2.boot[i, ] <- table(factor(ind.carriage.class.2.nonas[k], levels = c("carrier", "intermittent", "negative")))/sum(table(factor(ind.carriage.class.2.nonas[k], levels = c("carrier", "intermittent", "negative"))))
}

def1.carrier.cis <- quantile(carriage.1.boot[, 1], c(0.025, 0.975))
def1.intermittent.cis <- quantile(carriage.1.boot[, 2], c(0.025, 0.975))
def1.negative.cis <- quantile(carriage.1.boot[, 3], c(0.025, 0.975))

def2.carrier.cis <- quantile(carriage.2.boot[, 1], c(0.025, 0.975))
def2.intermittent.cis <- quantile(carriage.2.boot[, 2], c(0.025, 0.975))
def2.negative.cis <- quantile(carriage.2.boot[, 3], c(0.025, 0.975))

require(graphics)
# svg("./Plots/PropEwesInEachCarriageState_20160712.svg", height = 4, width = 3)
plot(x = 0, y = 0, xlim = c(.5, 3.5), ylim = c(0, 1), xaxt = "n",
     xlab = "", ylab = "Proportion of adult ewes", las = 1)
segments(x0 = .9, x1 = 1.1, y0 = carriage.status.props.1[1], y1 = carriage.status.props.1[1])
segments(x0 = 1.9, x1 = 2.1, y0 = carriage.status.props.1[2], y1 = carriage.status.props.1[2])
segments(x0 = 2.9, x1 = 3.1, y0 = carriage.status.props.1[3], y1 = carriage.status.props.1[3])
segments(x0 = 1, x1 = 1, y0 = def1.carrier.cis[1], y1 = def1.carrier.cis[2])
segments(x0 = 2, x1 = 2, y0 = def1.intermittent.cis[1], y1 = def1.intermittent.cis[2])
segments(x0 = 3, x1 = 3, y0 = def1.negative.cis[1], y1 = def1.negative.cis[2])
axis(side = 1, at = c(1, 2, 3), labels = c("Carrier", "Intermittent", "Negative"))

segments(x0 = 1.2, x1 = 1.4, y0 = carriage.status.props.2[1], y1 = carriage.status.props.2[1], col = "grey60")
segments(x0 = 2.2, x1 = 2.4, y0 = carriage.status.props.2[2], y1 = carriage.status.props.2[2], col = "grey60")
segments(x0 = 3.2, x1 = 3.4, y0 = carriage.status.props.2[3], y1 = carriage.status.props.2[3], col = "grey60")
segments(x0 = 1.3, x1 = 1.3, y0 = def2.carrier.cis[1], y1 = def2.carrier.cis[2], col = "grey60")
segments(x0 = 2.3, x1 = 2.3, y0 = def2.intermittent.cis[1], y1 = def2.intermittent.cis[2], col = "grey60")
segments(x0 = 3.3, x1 = 3.3, y0 = def2.negative.cis[1], y1 = def2.negative.cis[2], col = "grey60")
axis(side = 1, at = c(1, 2, 3), labels = c("Carrier", "Intermittent", "Negative"))
leg.text <- c("Definition 1: Positive in two consecutive years", "Definition 2: Positive twice in same year")
legend("top", leg.text, col = c("black", "grey60"), bty = "n", lty = c(1, 1))
# dev.off()