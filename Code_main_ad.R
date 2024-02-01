library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
source("Functions_ad.R")

init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5,2.6,1.8,2),eta=c(1.8,0.9,0.5,-2.2),sd=c(1.1,1.4,0.75,1, 0.5),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-1)
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho","theta_1","theta_2")

samsize= c(250, 500, 1000)

# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W
# the order below also matches the design order in section 4 (simulation study) of the paper

for(l in samsize)
{
  nsim = 100
  myseed = 750751
  message("sample size = ",l)
  #SimulationCI11_SaraIlias(l,nsim,myseed, init.value.theta_1, init.value.theta_2) # Design 1
  #SimulationCI12_SaraIlias(l,nsim,myseed, init.value.theta_1, init.value.theta_2) # Design 2
  SimulationCI21_SaraIlias(l,nsim,myseed, init.value.theta_1, init.value.theta_2) # Design 3
  SimulationCI22_SaraIlias(l,nsim,myseed,init.value.theta_1, init.value.theta_2) # Design 4
}


