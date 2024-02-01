
# Clear workspace
rm(list = ls())
# Load necessary packages, functions and data
fulldata <- read.csv("hiip_data_041515.csv")
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(mnormt)
library(matrixcalc)
library(cmprsk)
library(brglm2)
source("Functions_cmprsk.R")
source(paste(dirname(getwd()),"Functions_ad.R",sep="/"))

# data
dataset<-fulldata[(fulldata$followup_days != "M" & fulldata$followup_days != "I"),]
dataset$followup_days = as.numeric(dataset$followup_days)
dataset = dataset[dataset$followup_days>0,]
dataset <- dataset[1:50000,] # More observations results in memory error when fitting the model
dataset$age_entry = (dataset$age_entry-mean(dataset$age_entry))/sd(dataset$age_entry)
Y = as.matrix(log(dataset$followup_days))
Delta = as.matrix(dataset$final_vital_stat)
d1 = as.numeric(dataset$final_vital_stat==1) #death from breast cancer
d2 = as.numeric (dataset$final_vital_stat ==2) #death from other cause
da = as.numeric (dataset$final_vital_stat ==0) # alive or lost to follow up
dataset$intercept=rep(1,nrow(dataset))

Z = as.matrix(dataset$cohort==1) # selected + participated
W = as.matrix(dataset$cohort != 3) # selected
X = as.matrix(subset(dataset, select = c(ncol(dataset),age_entry)))
median(dataset$age_entry)

XandW = as.matrix(cbind(X,W))
n=nrow(dataset)
data=as.matrix(cbind(Y,d1,d2,X,Z,W))

namescoef =  c("beta_{T1,0}","beta_{T1,1}", "alpha_{T1}","lambda_{T1}",
               "beta_{T2,0}","beta_{T2,1}", "alpha_{T2}","lambda_{T2}",
               "sigma_{T1}","sigma_{T2}","rho","theta_1","theta_2")


# Define some useful variables
parl=ncol(X)+2
totparl=2*parl


# fitting the model
init.value.theta_1 <- 1
init.value.theta_2 <- 1


# Method using Firth
parlgamma=parl-1
results<-DataApplication_cmprsk.Firth(data, init.value.theta_1, init.value.theta_2) 

# Method using V=0 if W=0
parlgamma=parl-2
results<-DataApplication_cmprsk.alternative(data, init.value.theta_1, init.value.theta_2)

# results parametric estimation
results.screened <- results[[1]]
results.not_screened <- results[[2]]
results.control <- results[[3]]


# Nonparametric estimator of CIF
#datasub <- dataset[dataset$age_entry==50,]
datasub<-dataset

Y.sub = datasub$followup_days
d1.sub = as.numeric(datasub$final_vital_stat==1) #death from breast cancer
d2.sub = as.numeric (datasub$final_vital_stat ==2) #death from other cause
da.sub = as.numeric (datasub$final_vital_stat ==0) # alive or lost to follow up

D = d1.sub+2*d2.sub #1 of C1, 2 if C2 and 0 if censored (administrative)
Z.sub = as.matrix(datasub$cohort==1) # selected + participated
W.sub = as.matrix(datasub$cohort != 3) # selected
group=factor(Z.sub+2*W.sub) #0 (Z=0,W=0), 1(Z=1,W=0), 2(Z=0,W=1), 3(Z=1,W=1)

Time <- seq(from=1,to=8000,by=1)  # Times to evaluate
fit = timepoints(cuminc(Y.sub,D,group,cencode=0),Time) # nonparametric estimation

# Add these estimates to the estimates obtained from parametric estimation
results.screened <- rbind(results.screened, fit$est[3,], fit$est[6,])
results.not_screened <- rbind(results.not_screened, fit$est[2,], fit$est[5,])
results.control <- rbind(results.control, fit$est[1,], fit$est[4,])


# Plot results

par(mfrow=c(1,2))
# Death from breast cancer (C1)
plot(c(0,results.screened[1,]),c(0,results.screened[2,]),type="l", lty=3, main="Cancer",xlab="Time",ylab="Probability",xlim=c(0,8000), ylim=c(0,0.02)) # two-step estimator
lines(c(0,results.screened[1,]),c(0,results.screened[6,]), lty=3,type="l") # nonparametric estimator

lines(c(0,results.screened[1,]),c(0,results.not_screened[2,]), lty=1,type="l", col="grey60")
lines(c(0,results.screened[1,]),c(0,results.not_screened[6,]), lty=1,type="l", col="grey60")

lines(c(0,results.screened[1,]),c(0,results.control[2,]), lty=1,type="l")
lines(c(0,results.screened[1,]),c(0,results.control[6,]), lty=1,type="l")

# Death from other causes (C2)
plot(c(0,results.screened[1,]),c(0,results.screened[3,]),type="l", lty=3, main="Other cause",xlab="Time",ylab="Probability",xlim=c(0,8000), ylim=c(0,0.25))
lines(c(0,results.screened[1,]),c(0,results.screened[7,]), lty=3,type="l")

lines(c(0,results.screened[1,]),c(0,results.not_screened[3,]), lty=1,type="l", col="grey60")
lines(c(0,results.screened[1,]),c(0,results.not_screened[7,]), lty=1,type="l", col="grey60")

lines(c(0,results.screened[1,]),c(0,results.control[3,]), lty=1,type="l")
lines(c(0,results.screened[1,]),c(0,results.control[7,]), lty=1,type="l")


par(mfrow=c(1,1))

