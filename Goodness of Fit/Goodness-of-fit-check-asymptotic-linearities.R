
rm(list = ls())

library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(pbivnorm)
library(LaplacesDemon)
library(survival)
library(foreach)
library(doParallel)
library(expint)
library(fMultivar)

parent_directory <- dirname(getwd())
source(paste0(parent_directory, "/Functions_ad.R"))
source("Goodness-of-fit-test_functions.R")
source(paste0(parent_directory, "/Misspecification/Misspecification_functions.R"))

init.value.theta_1=1
init.value.theta_2=1
parN = list(beta=c(2.5, 2.6, 1.8, 2),
            eta=c(1.8, 0.9, 0.5, -2.2),
            sd=c(1.1, 1, 0.75, 1, 0.5),
            gamma=c(-1, 0.6, 2.3))
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = parl - 1
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho","theta_1","theta_2")


n <- 1000
iseed <- 5656
Zbin <- 1
Wbin <- 1
k_range_size <- 100
ltcm_estimation_size <- 5000
display.plot <- TRUE
nruns <- 10

data <- dat.sim.reg.debug(n, parN, iseed, Zbin, Wbin)
colnames(data) <- c("Y", "delta", "xi", "intercept", "x1", "Z", "W", "realV")

#### Estimate gamma and delta and useful vectors/matrices of derivatives ####

Y = data[,1]
Delta = data[,2]
Xi = data[,3]
X = data[,(5:(parl+1))]
Z = data[,parl+2]
W = data[,parl+3]
realV <- data[,parl+4]
XandW = cbind(data[,4],X,W)
n <- length(Y)

doDelta <- TRUE

if (Zbin == 1) {
  lin.mod <- lm(Z ~ X + W)
  gammaest <- lin.mod$coefficients
  V <- Z - XandW %*% gammaest
  
  summary(lin.mod)
} else {
  gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
}

M = cbind(data[,4:(2+parl)],V)

MnoV = data[,4:(3+parl)]

init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                 eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution


initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
initd[length(initd) - 2] <- 0
# structure of 'initd':
#
# [1:4] beta (4 params) = First 4 params of parhat1
# [5:8] eta (4 params) = Next 4 params of parhat1
# [9]   sigma1 = parhat1[9]
# [10]  sigma2 = parhat1[10]
# [11]  rho = 0
# [12]  theta_1 = parhat1[11]
# [13]  theta_2 = parhat1[12]

parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution

parhatG = c(parhat,as.vector(gammaest))
# - beta (4 params) = First 4 params of parhat
# - eta (4 params) = Next 4 params of parhat
# - sigma1 = parhat[9]
# - sigma2 = parhat[10]
# - rho = parhat[11]
# - theta_1 = parhat[12]
# - theta_2 = parhat[13]
# - gamma = (intercept, gamma_X, gamma_W)

prodvec = XandW[,1]

for (i in 1:parlgamma) {
  for (j in 2:parlgamma) {
    if (i<=j){
      prodvec <- cbind(prodvec, XandW[,i]*XandW[,j])
    }
  }
}

if (Zbin == 1) {
  sumsecder = c(rep(0,ncol(prodvec)))
  for (i in 1:length(sumsecder)) {
    sumsecder[i]= -sum(prodvec[,i])
  }
  
  # M-matrix: second derivative of m(W,Z,gamma)
  # WM here -> M in paper
  WM = sumsecder[1:parlgamma]
  for (i in 2:parlgamma) {
    newrow <- sumsecder[c(i,(i+2):(i+parlgamma))]
    WM <- rbind(WM,newrow) 
  }
  WM <- (2/n)*WM
  
  # Inverse of M-matrix
  # WMI here -> M^{-1} in paper
  WMI = ginv(WM)
  
  # mi here -> h_m(W, Z, \gamma^*) in paper
  mi = c()
  
  for(i in 1:n){
    ##########################################################################
    # For debugging: Use realV instead of V here
    
    # newrow <- 2*V[i]%*%XandW[i,]
    newrow <- 2*realV[i]%*%XandW[i,]
    
    ##########################################################################
    mi = rbind(mi,newrow)
  }
  mi <- 2*mi
  
  mi <- t(mi)
}

# psii here -> \Psi in paper
psii = -WMI%*%mi

if (doDelta) {
  
  # Hgamma here -> H in paper
  if (Zbin == 1) {
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  }
  
  # H here -> H_\delta in paper
  H = Hgamma[1:length(initd),1:length(initd)]
  
  # HI here -> Inverse of H_\delta in paper
  HI = ginv(H)
  
  # All derivatives with respect to a parameter of \delta AND \gamma.
  # Vargamma here -> H_\gamma in paper
  Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
  
  # gi here -> h_l(S_i, gamma, delta) in paper
  gi = c()
  
  for (i in 1:n)
  {
    J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    gi = rbind(gi,c(J1))
  }
  
  gi = t(gi)
  
  # Computation of \Sigma_\delta:
  
  # h_l(S, gamma, delta) + H_\gamma %*% Psi
  partvar = gi + Vargamma%*%psii
  
  # E[(h_l(S, gamma, delta) + H_\gamma %*% Psi)(h_l(S, gamma, delta) + H_\gamma %*% Psi)^T]
  Epartvar2 = (partvar%*%t(partvar))
  
  # H_\delta^{-1} E[...] (H_\delta^{-1})^T
  totvarex = HI%*%Epartvar2%*%t(HI)
}

#### Some useful variables, functions, etc. ####

Y_trans_theta1 <- YJtrans(Y, parhat[totparl+4])
Y_trans_theta2 <- YJtrans(Y, parhat[totparl+5])

# Asymptotic linear representation of \delta
Omega_delta <- - HI %*% (gi + Vargamma %*% psii)

# Asymptotic linear representation of \gamma
Omega_gamma <- psii

##############################################################################
# For debugging: Check if Omega_delta and Omega_gamma indeed do what they are
# supposed to do. I.e. satisfy lemma 2

gammaest
all(V[1:100] - (Z - XandW %*% gammaest)[1:100] < 10^(-10))

cbind(XandW %*% parN[["gamma"]] + V, Z)[1:100,]

# gamma_hat <- ginv(t(XandW) %*% XandW) %*% t(XandW) %*% (XandW %*% parN[["gamma"]] + realV)
gamma_hat <- ginv(t(XandW) %*% XandW) %*% t(XandW) %*% (Z)
gamma_hat - parN[["gamma"]]
ginv(t(XandW) %*% XandW) %*% t(XandW) %*% realV

influence_function <- matrix(rep(0, 3*n), ncol = 1000)
Sigma_X <- (1/n) * (t(XandW) %*% XandW)
inv_Sigma_X <- ginv(Sigma_X)
for (i in 1:n) {
  influence_function[,i] <- inv_Sigma_X %*% XandW[i,] %*% realV[i]
}

rowMeans(influence_function)

# \hat{\gamma} - \gamma^*
gamma_hat_minus_gamma <- gammaest - c(parN[["gamma"]])

data.frame("True diff" = gamma_hat_minus_gamma,
           "computed diff" = (1/n)*rowSums(Omega_gamma))

delta_hat_minus_delta <- parhat - c(parN[["beta"]], parN[["eta"]], parN[["sd"]])


data.frame("True diff" = delta_hat_minus_delta,
           "computed diff" = (1/n)*rowSums(Omega_delta))







