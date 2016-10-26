# Script to fit quadratic model of age by prevalence
samples <- read.csv("./Data/LostineSample2016Final_WithClasses_20160927.csv", header = T)
samples$movi.pcr.pos <- ifelse(samples$movi_qpcr == "Detected", 1,
                               ifelse(samples$movi_qpcr == "Not detected", 0, NA))

samples$animal.current.age <- samples$age_capture

for(i in 1:dim(samples)[1]){
samples$animal.current.age.2yrbins[i] <- ifelse(samples$animal.current.age[i] %in% c(0, 1), .5,
                                      ifelse(samples$animal.current.age[i] %in% c(2, 3), 2.5,
                                      ifelse(samples$animal.current.age[i] %in% c(4, 5), 4.5,
                                      ifelse(samples$animal.current.age[i] %in% c(6, 7), 6.5, 
                                      ifelse(samples$animal.current.age[i] %in% c(8, 9), 8.5,
                                      ifelse(samples$animal.current.age[i] %in% c(10, 11), 10.5, 
                                      ifelse(samples$animal.current.age[i] %in% c(12, 13), 12.5,
                                      ifelse(samples$animal.current.age[i] %in% c(14, 15), 14.5, 
                                      ifelse(samples$animal.current.age[i] %in% c(16, 17), 16.5,
                                      ifelse(samples$animal.current.age[i] %in% c(18, 19), 18.5,
                                      ifelse(samples$animal.current.age[i] %in% c(20), 20, NA)))))))))))
}

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