
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

n <- 5000
iseed <- 5656
Zbin <- 1
Wbin <- 1
k_range_size <- 100
ltcm_estimation_size <- 5000
display.plot <- TRUE
nruns <- 10
doParallel <- TRUE

#
# Some useful functions
#

gk <- function (x, Y, QK, dQ1k.bar) {
  # sum( ((1-QK)^(-2)*dQ1k.bar)[which(sort(Y) >= 0 & sort(Y) <= x)] )
  sum( ((1-QK)^(-2)*dQ1k.bar)[which(sort(Y) <= x)] )
}

dat.sim.reg.debug = function(n, par, iseed, Zbin, Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
  # bivariate normal distribution of error terms
  mu = c(0,0)
  sigma = matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
  err = mvrnorm(n, mu =mu , Sigma=sigma)
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2) {  # nu is standard logistic
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  } else if (Zbin==1) {# nu is standard normal
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T = IYJtrans(Mgen%*%beta+err1,sd[4]) # model time with real covariates
  C = IYJtrans(Mgen%*%eta+err2,sd[5]) # model censoring with real covariates
  
  par(mfrow = c(1, 2))
  min_T <- min(T)
  min_C <- min(C)
  hist(T, main = paste0("Minimum: ", min_T))
  hist(C, main = paste0("Minimum: ", min_C))
  
  A = runif(n,-20,20)
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  # nr of columns is nr of parameters
  # nr of rows is sample size
  
  Y = pmin(T,C,A) # observed non-transformed time
  d1 = as.numeric(Y==T) # censoring indicator
  xi1 = ifelse(Y==T,0,as.numeric(Y==C))
  data = cbind(Y,d1,xi1,M,realV) # data consisting of observed time,
  # censoring indicator, all data and the control function
  
  return(data)
}

get_surv_prob <- function(fit, times) {
  stepfun(fit$time, c(1, fit$surv))(times)
}

# Using the true parameters
GOF_debug <- function(data, k_range_size, ltcm_estimation_size) {
  
  #### Extract data ####
  
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  realV <- data[,parl+4]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  # NOTE: I changed 'V' here to 'realV'
  MrealV = cbind(data[,4:(2+parl)],realV)
  
  #### Some useful functions ####
  
  Y_trans_theta1_true <- YJtrans(Y, parN[["sd"]][4])
  Y_trans_theta2_true <- YJtrans(Y, parN[["sd"]][5])
  
  # Compute the cumulative density of K according to our two-step model
  FK_true <- rep(0, n)
  for (j in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1_true[j] - MrealV %*% parN[["beta"]])/parN[["sd"]][1]
    arg2_vector <- (Y_trans_theta2_true[j] - MrealV %*% parN[["eta"]])/parN[["sd"]][2]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK_true[j] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                        pbivnorm(cbind(arg1_vector, arg2_vector), rho = parN[["sd"]][3]))/n
  }
  FK_true <- FK_true[order(FK_true)]
  
  # Compute the censoring distribution of K. Since A ~ unif(0, 100), this will
  # simply be a the linear function through the origin and the point (100, 1)
  GK_true <- pmax(0, sort(Y)/40 + 0.5)
  
  par(mfrow = c(1, 1))
  
  plot(survfit(Surv(Y, (1 - Delta - Xi)) ~ 1, type = "kaplan-meier"))
  lines(sort(Y), 1 - GK_true, col = "red")
  
  QK_true <- 1 - (1 - FK_true) * (1 - GK_true)
  
  #### Define function to obtain evaluations of \Omega_KM####
  
  Omega_vct_KM <- function(k, dQ1k_true.bar, n) {
    
    Omegas_KM_true <- rep(0, n)
    
    part1 <- rep(0, n)
    part2 <- rep(0, n)
    part3 <- rep(0, n)
    
    for (i in 1:n) {
      
      FK.k_true <- stepfun(sort(Y), c(0, FK_true))(k)
      Omega_KM_true <- (1 - FK.k_true)*(gk(min(k, Y[i]), Y, QK_true, dQ1k_true.bar) + 
                                          ifelse((Y[i] <= k) & (Delta[i] + Xi[i] == 1),
                                                 1/(1 - QK_true[which(sort(Y) == Y[i])]),
                                                 0))
      # Omega_KM_true <- (gk(min(k, Y[i]), Y, QK_true, dQ1k_true.bar) + 
      #                                     ifelse((Y[i] <= k) & (Delta[i] + Xi[i] == 1),
      #                                            1/(1 - QK_true[which(sort(Y) == Y[i])]),
      #                                            0))
      
      
      Omegas_KM_true[i] <- Omega_KM_true
      
      part1[i] <- 1 - FK.k_true
      part2[i] <- gk(min(k, Y[i]), Y, QK_true, dQ1k_true.bar)
      part3[i] <- ifelse((Y[i] <= k) & (Delta[i] + Xi[i] == 1),
                         1/(1 - QK_true[which(sort(Y) == Y[i])]),
                         0)
    }
    
    list(Omegas_KM_true, part1, part2, part3)
  }
  
  #### Estimate variance of W_n(k) ####
  
  U <- min(max(Y[Delta == 1 | Xi == 1]), max(Y[Delta == 0 & Xi == 0]))
  
  K_RANGE_SIZE <- k_range_size
  
  # Determine range of k values in a data driven way
  k_range <- quantile(Y, seq(0, (K_RANGE_SIZE - 1))/(K_RANGE_SIZE - 1))
  
  if (doParallel) {
    if (!getDoParRegistered()) {
      message("No parallel backend registered. Continuing with sequential approach.")
      doParallel <- FALSE
    }
  }
  
  if (doParallel) {
    export.vector <- c("K_RANGE_SIZE", "parl", "totparl", "parlgamma")
    package.vector <- c("fMultivar")
    mean_var_estimates <- foreach(i = 1:K_RANGE_SIZE,
                                  .packages = package.vector,
                                  .export = export.vector,
                                  .combine = 'rbind') %dopar%
      {
        source("Goodness-of-fit-test_functions.R")
        k <- k_range[i]
        
        # Compute Q1(k) and dQ1
        Q1k_true <- rep(0, n)
        for (j in 1:n) {
          Q1k_true[j] <- sum( (Y < sort(Y)[j]) & (Delta + Xi == 1) )/n
        }
        Q1k_true.bar <- 1 - Q1k_true
        
        dQ1k_true.bar <- rep(0, n)
        dQ1k_true.bar[1] <- 1 - Q1k_true.bar[1]
        dQ1k_true.bar[2:n] <- Q1k_true.bar[2:n] - Q1k_true.bar[1:(n-1)]
        
        out <- Omega_vct_KM(k, dQ1k_true.bar, n)
        
        Omegas_KM_true <- out[[1]]
        part1 <- out[[2]]
        part2 <- out[[3]]
        part3 <- out[[4]]
        
        var.Omegas_KM_true1 <- var(Omegas_KM_true)
        var.Omegas_KM_true2 <- (1/n)*sum(Omegas_KM_true^2)
        
        var.part1 <- var(part1)
        var.part2 <- var(part2)
        var.part3 <- var(part3)
        
        cov.12 <- cov(part2, part3)
        
        c(mean(Omegas_KM_true), var.Omegas_KM_true1, var.part1, var.part2, var.part3,
          cov.12, var.Omegas_KM_true2, mean(part1), mean(part2), mean(part3))
      }
    
    mean_k <- mean_var_estimates[,1]
    var_k.KM_true <- mean_var_estimates[,2]
    
    var.part1_k <- mean_var_estimates[,3]
    var.part2_k <- mean_var_estimates[,4]
    var.part3_k <- mean_var_estimates[,5]
    
    cov.12_k <- mean_var_estimates[,6]
    var_k.KM_true2 <- mean_var_estimates[,7]
    
    mean.part1_k <- mean_var_estimates[,8]
    mean.part2_k <- mean_var_estimates[,9]
    mean.part3_k <- mean_var_estimates[,10]
  }
  
  #### Make some plots ####
  
  odd_points <- c(1, 2, 4, 5, 6, 8, 9, 13, 14, 17, 18, 19, 21, 22, 25, 27, 28,
                  30, 33, 37, 38, 43, 44, 46, 54, 55, 67, 81, 83)
  
  par(mfrow = c(1, 2))
  plot(k_range, mean_k)
  lines(k_range[odd_points], mean_k[odd_points], type = 'p', col = "red")
  plot(k_range, var_k.KM_true)
  
  
  
  par(mfrow = c(1, 2))
  plot(k_range, var_k.KM_true)
  plot(k_range, var_k.KM_true2)
  
  par(mfrow = c(2, 2))
  plot(k_range, mean_k)
  plot(k_range, mean.part1_k)
  plot(k_range, mean.part2_k)
  plot(k_range, mean.part3_k)
  
  
  # Plot the estimate of the st. dev. of the KM curve according to the survfit-
  # function and our theoretical result.
  par(mfrow = c(1, 1))
  KM <- survfit(Surv(Y, Delta + Xi) ~ 1, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain")
  plot(sort(Y), KM$std.err)
  
  # NOTE:
  # For some odd reason, std.err does not give the standard deviation based on
  # the greenwood formula but rather gives the same results as if you were to 
  # run this function with the commented \Omega_KM formula above (i.e. without
  # the (1 - FK.k_true) part in front)...
  
  lines(k_range, sqrt(var_k.KM_true)/sqrt(n), type = 'p', col = 'red')
  
  plot(KM)
  
  # Plot the variances of the different parts
  LAST_K <- 95
  
  par(mfrow = c(2, 2))
  plot(k_range[1:LAST_K], cov.12_k[1:LAST_K])
  plot(k_range[1:LAST_K], var.part1_k[1:LAST_K])
  plot(k_range[1:LAST_K], var.part2_k[1:LAST_K])
  plot(k_range[1:LAST_K], var.part3_k[1:LAST_K])
  
  # Greenwood's formula
  n_deaths <- rep(0, n)
  n_alive <- rep(0, n)
  
  for (i in 1:n) {
    n_alive[i] <- sum(Y >= sort(Y)[i])
    n_deaths[i] <- sum(Y == sort(Y)[i] & Delta + Xi == 1)
  }
  arg_sum <- n_deaths / (n_alive*(n_alive - n_deaths))
  
  greenwood_var <- rep(0, n)
  for (i in 1:n) {
    greenwood_var[i] <- get_surv_prob(KM, sort(Y)[i])^2 * sum(arg_sum[1:i])
  }
  
  par(mfrow = c(1, 1))
  plot(sort(Y), greenwood_var, ylim = c(0, max(c(greenwood_var, var_k.KM_true/n), na.rm = TRUE)))
  lines(k_range, var_k.KM_true/n, type = 'p', col = 'red')
}


# Using the estimated parameters
GOF_debug_2 <- function(data, k_range_size, ltcm_estimation_size) {

  #### Extract data and compute stuff ####
  
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  realV <- data[,parl+4]
  XandW = cbind(data[,4],X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    lin.mod <- lm(Z ~ X + W)
    gammaest <- lin.mod$coefficients
    V <- lin.mod$residuals
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
  
  #### Some useful functions ####
  
  Y_trans_theta1 <- YJtrans(Y, parhat[totparl+4])
  Y_trans_theta2 <- YJtrans(Y, parhat[totparl+5])
  
  # Compute the cumulative density of K according to our two-step model
  FK <- rep(0, n)
  for (j in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[j] - M %*% parhat[1:parl])/parhat[totparl+1]
    arg2_vector <- (Y_trans_theta2[j] - M %*% parhat[(parl+1):totparl])/parhat[totparl+2]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[j] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[totparl+3]))/n
  }
  FK <- FK[order(FK)]
  
  # Compute the censoring distribution of K
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1, type = "kaplan-meier")
  GK <- 1 - get_surv_prob(surv_A, sort(Y))
  
  # Compute QK (useful when computing \Omega_{KM})
  QK <- 1 - (1 - FK) * (1 - GK)
  
  #### Define function to obtain evaluations of \Omega_KM####
  
  Omega_vct_KM <- function(k, dQ1k.bar, n) {
    
    Omegas_KM <- rep(0, n)
    
    part1 <- rep(0, n)
    part2 <- rep(0, n)
    part3 <- rep(0, n)
    
    for (i in 1:n) {
      
      FK.k <- stepfun(sort(Y), c(0, FK))(k)
      Omega_KM <- (1 - FK.k)*(gk(min(k, Y[i]), Y, QK, dQ1k.bar) + 
                                          ifelse((Y[i] <= k) & (Delta[i] + Xi[i] == 1),
                                                 1/(1 - QK[which(sort(Y) == Y[i])]),
                                                 0))
      # Omega_KM_true <- (gk(min(k, Y[i]), Y, QK_true, dQ1k_true.bar) + 
      #                                     ifelse((Y[i] <= k) & (Delta[i] + Xi[i] == 1),
      #                                            1/(1 - QK_true[which(sort(Y) == Y[i])]),
      #                                            0))
      
      
      Omegas_KM[i] <- Omega_KM
      
      part1[i] <- 1 - FK.k
      part2[i] <- gk(min(k, Y[i]), Y, QK, dQ1k.bar)
      part3[i] <- ifelse((Y[i] <= k) & (Delta[i] + Xi[i] == 1),
                         1/(1 - QK[which(sort(Y) == Y[i])]),
                         0)
    }
    
    list(Omegas_KM, part1, part2, part3)
  }
  
  #### Estimate variance of W_n(k) ####
  
  U <- min(max(Y[Delta == 1 | Xi == 1]), max(Y[Delta == 0 & Xi == 0]))
  
  K_RANGE_SIZE <- k_range_size
  
  # Determine range of k values in a data driven way
  k_range <- quantile(Y, seq(0, (K_RANGE_SIZE - 1))/(K_RANGE_SIZE - 1))
  
  if (doParallel) {
    if (!getDoParRegistered()) {
      message("No parallel backend registered. Continuing with sequential approach.")
      doParallel <- FALSE
    }
  }
  
  if (doParallel) {
    export.vector <- c("K_RANGE_SIZE", "parl", "totparl", "parlgamma")
    package.vector <- c("fMultivar")
    mean_var_estimates <- foreach(i = 1:K_RANGE_SIZE,
                                  .packages = package.vector,
                                  .export = export.vector,
                                  .combine = 'rbind') %dopar%
      {
        source("Goodness-of-fit-test_functions.R")
        k <- k_range[i]
        
        # Compute Q1(k) and dQ1
        Q1k <- rep(0, n)
        for (j in 1:n) {
          Q1k[j] <- sum( (Y < sort(Y)[j]) & (Delta + Xi == 1) )/n
        }
        Q1k.bar <- 1 - Q1k
        
        dQ1k.bar <- rep(0, n)
        dQ1k.bar[1] <- 1 - Q1k.bar[1]
        dQ1k.bar[2:n] <- Q1k.bar[2:n] - Q1k.bar[1:(n-1)]
        
        out <- Omega_vct_KM(k, dQ1k.bar, n)
        
        Omegas_KM <- out[[1]]
        part1 <- out[[2]]
        part2 <- out[[3]]
        part3 <- out[[4]]
        
        var.Omegas_KM <- var(Omegas_KM)
        
        var.part1 <- var(part1)
        var.part2 <- var(part2)
        var.part3 <- var(part3)
        
        cov.12 <- cov(part2, part3)
        
        c(mean(Omegas_KM), var.Omegas_KM, var.part1, var.part2, var.part3,
          cov.12, mean(part1), mean(part2), mean(part3))
      }
    
    mean_k <- mean_var_estimates[,1]
    var_k.KM <- mean_var_estimates[,2]
    
    var.part1_k <- mean_var_estimates[,3]
    var.part2_k <- mean_var_estimates[,4]
    var.part3_k <- mean_var_estimates[,5]
    
    cov.12_k <- mean_var_estimates[,6]
    
    mean.part1_k <- mean_var_estimates[,7]
    mean.part2_k <- mean_var_estimates[,8]
    mean.part3_k <- mean_var_estimates[,9]
  }
  
  #### Make some plots ####
  
  odd_points <- c(1, 2, 4, 5, 6, 8, 9, 13, 14, 17, 18, 19, 21, 22, 25, 27, 28,
                  30, 33, 37, 38, 43, 44, 46, 54, 55, 67, 81, 83)
  
  par(mfrow = c(1, 2))
  plot(k_range, mean_k)
  lines(k_range[odd_points], mean_k[odd_points], type = 'p', col = "red")
  plot(k_range, var_k.KM)
  
  par(mfrow = c(2, 2))
  plot(k_range, mean_k)
  plot(k_range, mean.part1_k)
  plot(k_range, mean.part2_k)
  plot(k_range, mean.part3_k)
  
  
  # Plot the estimate of the st. dev. of the KM curve according to the survfit-
  # function and our theoretical result.
  par(mfrow = c(1, 1))
  KM <- survfit(Surv(Y, Delta + Xi) ~ 1, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain")
  plot(sort(Y), KM$std.err)
  lines(k_range, sqrt(var_k.KM)/sqrt(n), type = 'p', col = 'red')
  
  plot(KM)
  
  # Plot the variances of the different parts
  LAST_K <- 95
  
  par(mfrow = c(2, 2))
  plot(k_range[1:LAST_K], cov.12_k[1:LAST_K])
  plot(k_range[1:LAST_K], var.part1_k[1:LAST_K])
  plot(k_range[1:LAST_K], var.part2_k[1:LAST_K])
  plot(k_range[1:LAST_K], var.part3_k[1:LAST_K])
  
  # Greenwood's formula
  n_deaths <- rep(0, n)
  n_alive <- rep(0, n)
  
  for (i in 1:n) {
    n_alive[i] <- sum(Y >= sort(Y)[i])
    n_deaths[i] <- sum(Y == sort(Y)[i] & Delta + Xi == 1)
  }
  arg_sum <- n_deaths / (n_alive*(n_alive - n_deaths))
  
  greenwood_var <- rep(0, n)
  for (i in 1:n) {
    greenwood_var[i] <- get_surv_prob(KM, sort(Y)[i])^2 * sum(arg_sum[1:i])
  }
  
  par(mfrow = c(1, 1))
  plot(sort(Y), greenwood_var, ylim = c(0, max(c(greenwood_var, var_k.KM/n), na.rm = TRUE)))
  lines(k_range, var_k.KM/n, type = 'p', col = 'red')
}

#
# Generate data
#

data <- dat.sim.reg.debug(n, parN, iseed, Zbin, Wbin)
colnames(data) <- c("Y", "delta", "xi", "intercept", "x1", "Z", "W", "realV")

clust <- makeCluster(10)
registerDoParallel(clust)

#
# Main function
#

GOF_debug(data, k_range_size, ltcm_estimation_size)
GOF_debug_2(data, k_range_size, ltcm_estimation_size)


# When all times are positive and A ~ Unif(0, 100):

# It seems like comparing the results of survfits $std.err with the obtained
# results from Lo and Singh was just a red haring. Computing the Greenwood
# formula-variance ourselves we do obtain the correct result.

# NOTE: to make the times positive, we added 30 and 10 to the intercepts of the
# model for T and C, respectively.


# When T and C are as specified in all of our main code and A ~ Unif(-20, 20):

# Results still look very good. It can be seen, however, that a larger sample
# size is needed to get a good overlap of the greenwood variance and the Lo and
# Singh variance (n ~= 5000).

# When T and C are as specified in all of our main code, A ~ Unif(-20, 20) and
# we work on the estimated parameters instead of the true ones, we still get 
# good overlap of both curved when n = 5000.





