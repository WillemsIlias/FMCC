rm(list = ls())

library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(mnormt)
library(matrixcalc)
library(cmprsk)
source("Functions_cmprsk.R")

init.value.theta_1=1
init.value.theta_2=1
init.value.theta_3=1
parN = list(beta1=c(4,1,1.8,2),beta2=c(4.65,-1,-1.2,1.3),beta3=c(4.15,1.5,-1.5,-1.8),sd=c(1.5,1,1.2,-0.80,0.65,-0.75,1.5,0.5, 1),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 3*parl
parlgamma = (parl-1)
namescoef =  c("beta_{01}","beta_{11}","alpha_1","lambda_1","beta_{02}","beta_{12}","alpha_2","lambda_2","beta_{03}","beta_{13}","alpha_3","lambda_3","sigma_1","sigma_2","sigma_3","rho_12","rho_13","rho_23","theta_1","theta_2","theta_3")

samsize= c(1000)


# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W

nsim = 500
myseed = 750751

for(l in samsize)
{
  message("sample size = ",l)
  results=SimulationCI22_Sara(l,nsim,myseed,init.value.theta_1, init.value.theta_2, init.value.theta_3) 
}

########### Results of estimation

results.mean = results[[1]]
results.RMSE = results[[2]]


############## Plot curves

par(mfrow=c(2,2))
#C1-Z1-W1
plot(c(0,results.mean[1,]),c(0,results.mean[2,]),type="l", lty=3, main=expression(paste(C1-Z1-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,30), ylim=c(0,0.35)) # two-step model
lines(c(0,results.mean[1,]),c(0,results.mean[10,]), lty=3,type="l",col="grey50") # nonparametric
lines(c(0,results.mean[1,]),c(0,results.mean[6,]), lty=1,type="l", col="grey70") # naive
lines(c(0,results.mean[1,]),c(0,results.mean[8,]), lty=1,type="l") # true

#C2-Z1-W1
plot(c(0,results.mean[1,]),c(0,results.mean[3,]),type="l", lty=3, main=expression(paste(C2-Z1-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,30))
lines(c(0,results.mean[1,]),c(0,results.mean[11,]), lty=3,type="l",col="grey50")
lines(c(0,results.mean[1,]),c(0,results.mean[7,]), lty=1,type="l", col="grey70")
lines(c(0,results.mean[1,]),c(0,results.mean[9,]), lty=1,type="l")

#C1-Z0-W0
plot(c(0,results.mean[1,]),c(0,results.mean[32,]),type="l", lty=3, main=expression(paste(C1-Z0-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,100), ylim=c(0,0.8))
lines(c(0,results.mean[1,]),c(0,results.mean[40,]), lty=3,type="l",col="grey50")
lines(c(0,results.mean[1,]),c(0,results.mean[36,]), lty=1,type="l", col="grey70")
lines(c(0,results.mean[1,]),c(0,results.mean[38,]), lty=1,type="l")

#C2-Z0-W0
plot(c(0,results.mean[1,]),c(0,results.mean[33,]),type="l", lty=3, main=expression(paste(C2-Z0-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,100), ylim=c(0,0.05))
lines(c(0,results.mean[1,]),c(0,results.mean[41,]), lty=3,type="l",col="grey50")
lines(c(0,results.mean[1,]),c(0,results.mean[37,]), lty=1,type="l", col="grey70")
lines(c(0,results.mean[1,]),c(0,results.mean[39,]), lty=1,type="l")
par(mfrow=c(1,1))




par(mfrow=c(2,2))
#C1-Z0-W1
plot(c(0,results.mean[1,]),c(0,results.mean[22,]),type="l", lty=1, main=expression(paste(C1-Z0-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[24,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[26,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[28,]), lty=1,type="l", col="red")

#C2-Z0-W1
plot(results.mean[1,],results.mean[23,],type="l", lty=1, main=expression(paste(C2-Z0-tilde(W),1)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[25,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[27,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[29,]), lty=1,type="l", col="red")

#C1-Z1-W0
plot(c(0,results.mean[1,]),c(0,results.mean[12,]),type="l", lty=1, main=expression(paste(C1-Z1-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[14,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[16,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[18,]), lty=1,type="l", col="red")

#C2-Z1-W0
plot(c(0,results.mean[1,]),c(0,results.mean[13,]),type="l", lty=1, main=expression(paste(C2-Z1-tilde(W),0)), xlab="Time",ylab="Probability",xlim=c(0,50))
lines(c(0,results.mean[1,]),c(0,results.mean[15,]), lty=5,type="l")
lines(c(0,results.mean[1,]),c(0,results.mean[17,]), lty=4,type="l", col="azure4")
lines(c(0,results.mean[1,]),c(0,results.mean[19,]), lty=1,type="l", col="red")
par(mfrow=c(1,1))

# look at results
results.mean[,c(19,39,59,99,199)]
results.RMSE[,c(19,39,59,99,199)]


##################### Global RMSE
# approximate the global RMSE integral on a grid

# Define the grid
grid.lower.bound <- 1
grid.upper.bound1 <- 30
grid.upper.bound2 <- 100

number.of.grid.cells1 <- (grid.upper.bound1-grid.lower.bound)*2
number.of.grid.cells2 <- (grid.upper.bound2-grid.lower.bound)*2
cell.width <- 0.5

# Approximate the integrals for Z=1, W=1 
int.C1.11 <- 0
int.C2.11 <- 0
int.C1E.11 <- 0
int.C2E.11 <- 0
int.C1C.11 <- 0
int.C2C.11 <- 0

# RMSE 
for (i in 1:number.of.grid.cells1) {
  int.C1.11 <- int.C1.11 + results.RMSE[2,i]*(cell.width)
  int.C2.11 <- int.C2.11 + results.RMSE[3,i]*(cell.width)
  int.C1E.11 <- int.C1E.11 + results.RMSE[6,i]*(cell.width)
  int.C2E.11 <- int.C2E.11 + results.RMSE[7,i]*(cell.width)
  int.C1C.11 <- int.C1C.11 + results.RMSE[8,i]*(cell.width)
  int.C2C.11 <- int.C2C.11 + results.RMSE[9,i]*(cell.width)
}

# Approximate the integrals for Z=0, W=0
int.C1.00 <- 0
int.C2.00 <- 0
int.C1E.00 <- 0
int.C2E.00 <- 0
int.C1C.00 <- 0
int.C2C.00 <- 0

# RMSE 
for (i in 1:number.of.grid.cells2) {
  int.C1.00 <- int.C1.00 + results.RMSE[26,i]*(cell.width)
  int.C2.00 <- int.C2.00 + results.RMSE[27,i]*(cell.width)
  int.C1E.00 <- int.C1E.00 + results.RMSE[30,i]*(cell.width)
  int.C2E.00 <- int.C2E.00 + results.RMSE[31,i]*(cell.width)
  int.C1C.00 <- int.C1C.00 + results.RMSE[32,i]*(cell.width)
  int.C2C.00 <- int.C2C.00 + results.RMSE[33,i]*(cell.width)
}






