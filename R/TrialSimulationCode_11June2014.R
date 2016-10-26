#------------------------------------------------------#
#-- Code for analysis of sheep longitudinal sampling --#
#------------------------------------------------------#

#----------------------#
#-- 1) Simulate data --#
#----------------------#

#-- 1.1 Specify sample size 
N <- 25

#-- 1.2 Specify success prob (based on state) for each individual
p.success <- c(rep(.9, floor(N/3)), rep(0, N / 3), rep(.5, N / 3), .5)

#-- 1.3 Generate number of trials per individual
n.trials <- floor(runif(N, 3, 7))

#-- 1.4 Generate test outcomes for each individual
n.success <- rep(NA, N)
for(i in 1:N){
  n.success[i] <- rbinom(1, n.trials[i], prob = p.success[i])
}

#-- 1.5 Calculate distances between each binomial observation 
require(vegan)
#-- use Mountford or Raup-Crick indices for pres/abs data in vegdist function
dist.mat <- vegdist()

