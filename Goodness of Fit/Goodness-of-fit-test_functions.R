
################################################################################
#                             Some useful functions                            #
################################################################################

# to obtain quantile from survival function
get_surv_prob <- function(fit, times) {
  stepfun(fit$time, c(1, fit$surv))(times)
}

# Weight function for test statistic
weight_fun <- function(x, U) {
  # as.numeric((x > -U) & (x < U))
  as.numeric(x < U)
}

DthetaYJtrans <- function (y, theta) {
  if (theta == 0) {
    warning(paste("Given value for theta lies on boundary of parameter space.",
                  "Using theta = 0.00001 instead"))
    theta <- 0.00001
  }
  if (theta == 2) {
    warning(paste("Given value for theta lies on boundary of parameter space.",
                  "Using theta = 1.99999 instead"))
    theta <- 1.99999
  }
  
  DthetaYJtrans_elementwise <- function (elem, theta) {
    if (elem >= 0) {
      result <- ((elem+1)^theta * (theta*log(elem+1) - 1) + 1) / (theta^2)
    } else {
      result <- ((-elem+1)^(2 - theta) * ((2 - theta)*log(-elem+1) - 1) + 1) / ((2 - theta)^2)
    }
  }
  
  unlist(lapply(y, DthetaYJtrans_elementwise, theta = theta))
}

gk <- function (x, Y, QK, dQ1k.bar) {
  sum( ((1-QK)^(-2)*dQ1k.bar)[which(sort(Y) < x)] )
}

dat.sim.reg.debug = function(n,par,iseed,Zbin,Wbin){
  
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
  A = runif(n,-20,20)
  
  par(mfrow = c(1, 3))
  min_T <- min(T)
  min_C <- min(C)
  min_A <- min(A)
  hist(T, main = paste0("Minimum: ", min_T))
  hist(C, main = paste0("Minimum: ", min_C))
  hist(A, main = paste0("Minimum: ", min_A))
  
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

################################################################################
#                            Main GOF test functions                           #
################################################################################

#
# Functions bootstrap method
#

GOF_test_parallel <- function(data, B, iseed, Zbin, Wbin, display.plot,
                              multiple.starting.points = FALSE) {
  
  #### Input validation ####
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid specification of Zbin and/or Wbin (can only be 1 or 2)")
  }
  
  #### Compute test statistic ####
  
  # Estimate the parameter vectors gamma and delta
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  intercept = data[,4]
  X = data[,(5:(parl+1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  XandW = cbind(intercept,X,W)
  n <- length(Y)
  
  if (Zbin == 1) {
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
  # Maximization of the likelihood
  if (multiple.starting.points) {
    
    # If multiple starting points should be used for the gradient descent
    # algorithm...
    
    # Define starting values of the parameter vector
    coef.starts <- c(-3, 0, 3)
    theta.starts <- c(0.5, 1.5)
    
    par.starts <- list()
    for (coef.idx in 1:totparl) {
      for (coef.start in coef.starts) {
        for (theta1.start in theta.starts) {
          for (theta2.start in theta.starts) {
            par.starts[[length(par.starts) + 1]] <-
              c(coef.idx, coef.start, theta1.start, theta2.start)
          }
        }
      }
    }
    
    subset <- list(par.starts[[5]], par.starts[[6]], par.starts[[7]],
                   par.starts[[8]])
    par.starts <- subset
    
    # Initialize an object that will store all the parameter vectors
    parhats.list <- list()
    
    # For each starting vector, run the optimization
    for (par.start in par.starts) {
      
      # Define the initial parameter vector of this iteration
      init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
      init[par.start[1]] <- par.start[2]
      init[length(init) - 1] <- par.start[3]
      init[length(init)] <- par.start[4]
      
      # Run the optimization
      parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
      
      initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
      initd[length(initd) - 2] <- 0
      
      parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
      
      # Store the result
      parhats.list[[length(parhats.list) + 1]] <- list(
        init = init,
        parhat = parhat
      )
    }
    
    # Find the parameter vector corresponding to the pseudo-global maximum of
    # the likelihood function
    # (note that LikF returns the negative of the log-likelihood)
    min.lik <- Inf
    parhat <- NULL
    for (entry in parhats.list) {
      lik.eval <- LikF(entry[[2]], Y, Delta, Xi, M)
      if (min.lik > lik.eval) {
        min.lik <- lik.eval
        parhat <- entry[[2]]
      }
    }
    
  } else {
    
    # If only one starting point should be used for the gradient descent
    # algorithm...
    
    init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
  }
  
  # Compute F_K(k; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[totparl+4])
  Y_trans_theta2 <- YJtrans(Y, parhat[totparl+5])
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FK <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:parl])/parhat[totparl+1]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[(parl+1):totparl])/parhat[totparl+2]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[totparl+3]))/n
  }
  
  # Order the values
  FK <- FK[order(FK)]
  
  # Useful when computing the test statistic later on. In a sense, it can be
  # seen to replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function of K = min(T, C). Since
  # the censoring variable A is independent of K, this can simply be done using
  # a Kaplan-Meier estimator.
  K_delta <- as.numeric((Delta == 1) | (Xi == 1)) # Censoring indicator for K = min(T,C)
  surv_km <- survfit(Surv(Y, K_delta) ~ 1, type = "kaplan-meier")
  cdf_km <- 1 - get_surv_prob(surv_km, sort(Y))
  
  # If necessary, display the results
  if (display.plot) {
    old.mfrow.setting <- par()$mfrow
    par(mfrow = c(1, 2))
    plot(sort(Y), FK, type = 'l', main = "",
    xlab = "K = min(T,C)", ylab = "F(K)")
    lines(sort(Y), cdf_km, type = 's', col = "red")
    legend("topleft", legend = c("Estimated model", "Kaplan-Meier"),
           col = c("black", "red"), lty = 1)
  }
  
  # Compute the appropriate test statistic
  TCM <- n*sum((cdf_km - FK)^2 * w)
  
  #### Compute the bootstrap distribution ####
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1, type = "kaplan-meier")
  cdf_A = 1 - get_surv_prob(surv_A, Y[order(Y)])
  
  package.vector <- c("MASS", "nloptr", "VGAM", "pbivnorm", "survival")
  export.vector <- c("parl", "totparl", "parent_directory")
  TCMb_vector <- foreach(b = 1:B,
                         .packages = package.vector,
                         .export = export.vector,
                         .combine = c) %dopar% {

    source(paste0(parent_directory, "/Functions_ad.R"))
    source(paste0(parent_directory, "/Goodness of Fit/Goodness-of-fit-test_functions.R"))
    
    MAX_SEED_SIZE <- 2147483647
    
    # Create bootstrap data sample
    bootstrapseed <- (B*iseed + b) %% MAX_SEED_SIZE
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star <- rep(0, n)
    
    # If there is no administrative censoring in the data, then the Kaplan - 
    # Meier estimator 'surv_A' will always be 1. In this case, we set all
    # simulated administrative censoring times equal to infinity. Otherwise,
    # proceed normally.
    if (all(surv_A$surv == 1)) {
      A.star <- rep(Inf, n)
    } else {
      A.star <- quantile(surv_A, probs = runif(n, 1e-10, 1 - 1e-10))$quantile
      A.star <- as.numeric(A.star)
      A.star[is.na(A.star)] <- Inf
    }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star == T.star)
    Xi.star = ifelse(Y.star == T.star, 0, as.numeric(Y.star == C.star))
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[totparl+4])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[totparl+5])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FK.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:parl])/parhat.star[totparl+1]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[(parl+1):totparl])/parhat.star[totparl+2]
      
      # Evaluate them in standard normal CDF, and compute the average
      FK.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[totparl+3]))/n
    }
    
    FK.star <- FK.star[order(FK.star)]
    
    w.star = rep(0,n)
    w.star[1] = FK.star[1]
    w.star[2:n] = FK.star[2:n]-FK.star[1:(n-1)]
    
    K_delta.star <- as.numeric((Delta.star == 1) | (Xi.star == 1))
    surv_km.star <- survfit(Surv(Y.star, K_delta.star) ~ 1, type = "kaplan-meier")
    cdf_km.star <- 1 - get_surv_prob(surv_km.star, sort(Y.star))
    
    # Output of loop iteration
    n*sum((cdf_km.star - FK.star)^2 * w.star)
  } 
  
  #### Extract, plot and return results ####
  
  if (display.plot) {
    hist(TCMb_vector, main = "",
         xlab = "Bootstrapped test statistics",
         xlim = c(0, max(max(TCMb_vector, na.rm = TRUE) + 0.1, TCM)))
    abline(v = TCM, col = "red")
    
    # Reset old mfrow setting
    par(mfrow = old.mfrow.setting)
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9, na.rm = TRUE))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95, na.rm = TRUE))
  significant3 <- (TCM > quantile(TCMb_vector, prob=0.80, na.rm = TRUE))
  pvalue <- sum(TCM < TCMb_vector, na.rm = TRUE)/length(TCMb_vector)
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2, signif80 = significant3, pvalue = pvalue)
}


GOF_EstimateRunDuration <- function(data, B, nruns, Zbin, Wbin, parallel) {
  testseed <- 123321
  start.time <- Sys.time()
  if (!parallel) {
    out <- GOF_test(data, B, testseed, Zbin, Wbin, FALSE)
  } else {
    clust <- makeCluster(10)
    registerDoParallel(clust)
    
    out <- GOF_test_parallel(data, B, testseed, Zbin, Wbin, FALSE)
    
    stopCluster(clust)
  }
  
  diff <- difftime(Sys.time(), start.time, units = "secs")
  diff <- as.integer(round(diff))
  duration <- round((diff*nruns)/60)
  if (duration > 60) {
    message("The simulation will take approximately ", floor(duration / 60),
            " hours and ", duration %% 60, " minutes")
  } else {
    message("The simulation will take approximately ", duration, " minutes")
  } 
}


GOF_SimulationTypeIerror <- function(parN, nruns, B, iseed, Zbin, Wbin, parallel,
                                     save.bootstrap.distributions) {
  nbr_reject90 <- 0
  nbr_reject95 <- 0
  TCM_reject <- c()
  TCMb_reject <- c()
  
  # Set up for parallel computing if necessary
  if (parallel) {
    clust <- makeCluster(10)
    registerDoParallel(clust)
  }
  
  for (run in 1:nruns) {
    message(paste0("Currently busy with run ", run, " out of ", nruns, "."))
    run_seed <- iseed + run
    
    data <- dat.sim.reg(n, parN, run_seed, Zbin, Wbin)
    
    if (parallel) {
      out <- GOF_test_parallel(data, B, run_seed, Zbin, Wbin, display.plot = FALSE)
    } else {
      out <- GOF_test(data, B, run_seed, Zbin, Wbin, display.plot = FALSE)
    }
    
    
    if (out[[3]]) {
      nbr_reject90 <- nbr_reject90 + 1
      TCM_reject <- cbind(TCM_reject, out[[1]])
      TCMb_reject <- cbind(TCMb_reject, out[[2]])
    }
    
    if (out[[4]]) {
      nbr_reject95 <- nbr_reject95 + 1
    }
    
    if (save.bootstrap.distributions) {
      filename <- paste0(run_seed, ".csv")
      folder <- paste0("TypeIerror/BootstrapTCM_B", B, "_n", n)
      
      # If folder doesn't exist, create it.
      if (!file.exists(folder)) {
        dir.create(folder)
      }
      
      write.csv(out[[2]], file = paste0(folder, "/", filename))
    }
  }
  
  # stop cluster after parallel computing
  if (parallel) {
    stopCluster(clust)
  }
  
  list(reject90 = nbr_reject90/nruns, reject95 = nbr_reject95/nruns,
       TCM_rejected = TCM_reject, TCMb_rejected = TCMb_reject)
}


GOF_SimulationMisspecification <- function(type, par, nruns, B, iseed, Zbin, Wbin, parallel) {
  nbr_reject90 <- 0
  nbr_reject95 <- 0
  TCM_reject <- c()
  TCMb_reject <- c()
  
  # Set up for parallel computing if necessary
  if (parallel) {
    clust <- makeCluster(10)
    registerDoParallel(clust)
  }
  
  for (run in 1:nruns) {
    run_seed <- iseed + run
    message(paste0("Currently busy with run ", run, " out of ", nruns, "."))
    
    if (parallel) {
      if (type == "skew") {
        data.skew <- data.misspecified.skew(n, par, run_seed, Zbin, Wbin)[[1]]
        out <- GOF_test_parallel(data.skew, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "t") {
        data.t <- data.misspecified.t(n, par, run_seed, Zbin, Wbin)
        out <- GOF_test_parallel(data.t, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "heteroscedastic") {
        data.heteroscedastic <- data.misspecified.heteroscedastic(n, par, run_seed, Zbin, Wbin)
        out <- GOF_test_parallel(data.heteroscedastic, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'probit') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as probit control function requires ", 
                         "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.probit <- data.misspecified.probit(n, par, run_seed, Wbin)
        out <- GOF_test_parallel(data.probit, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'cloglog') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as cloglog control function requires ", 
                         "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.probit <- data.misspecified.cloglog(n, par, run_seed, Wbin)
        out <- GOF_test_parallel(data.probit, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else {
        stop(paste0("'type' parameter must be either 'skew', 't', 'heteroscedastic'",
                    ", 'probit' or 'cloglog'."))
      }
    } else {
      if (type == "skew") {
        data.skew <- data.misspecified.skew(n, par, run_seed, Zbin, Wbin)[[1]]
        out <- GOF_test(data.skew, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "t") {
        data.t <- data.misspecified.t(n, par, run_seed, Zbin, Wbin)
        out <- GOF_test(data.t, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == "heteroscedastic") {
        data.heteroscedastic <- data.misspecified.heteroscedastic(n, par, run_seed, Zbin, Wbin)
        out <- GOF_test(data.heteroscedastic, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'probit') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as probit control function requires ", 
                         "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.probit <- data.misspecified.probit(n, par, run_seed, Wbin)
        out <- GOF_test(data.probit, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else if (type == 'cloglog') {
        if (Zbin == 1) {
          warning(paste0("Argument Zbin = 1 invalid as cloglog control function requires ", 
                         "binary endogenous variable. Using Zbin = 2 instead."))
          Zbin <- 2
        }
        data.cloglog <- data.misspecified.cloglog(n, par, run_seed, Wbin)
        out <- GOF_test(data.cloglog, B, run_seed, Zbin, Wbin, display.plot = FALSE)
        
      } else {
        stop(paste0("'type' parameter must be either 'skew', 't', 'heteroscedastic'",
                    ", 'probit' or 'cloglog'."))
        
      }
    }
    
    if (out[[3]]) {
      nbr_reject90 <- nbr_reject90 + 1
      TCM_reject <- cbind(TCM_reject, out[[1]])
    }
    if (out[[4]]) {
      nbr_reject95 <- nbr_reject95 + 1
      TCMb_reject <- cbind(TCMb_reject, out[[2]])
    }
  }
  
  # stop cluster after parallel computing
  if (parallel) {
    stopCluster(clust)
  }
  
  list(reject90 = nbr_reject90/nruns, reject95 = nbr_reject95/nruns,
       TCM_rejected = TCM_reject, TCMb_rejected = TCMb_reject)
}


#
# Functions asymptotic theory
#

GOF_test_noBootstrap <- function(data, Zbin, Wbin, k_range_size, ltcm_estimation_size, 
                                 display.plot) {
  
  # Since this method still doesn't work as it should, throw an error an error
  # every time it is used. We leave the code below so that it can be fixed at a
  # later moment in time.
  
  stop(paste0("The goodness-of-fit test based on the limiting distribution is ",
              "not (yet) correctly implemented. Use the bootstrap version ",
              "instead"))
  
  verbose <- TRUE # Maybe add this as extra functionality later
  doParallel <- TRUE # Maybe add this as extra functionality later
  
  if ((Zbin != 1 && Zbin != 2) | (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  if (Zbin == 2) {
    stop("Not implemented yet for Zbin = 2")
  }
  
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
  # - beta (4 params) = First 4 params of parhat
  # - eta (4 params) = Next 4 params of parhat
  # - sigma1 = parhat[9]
  # - sigma2 = parhat[10]
  # - rho = parhat[11]
  # - theta_1 = parhat[12]
  # - theta_2 = parhat[13]
  # - gamma = (intercept, gamma_X, gamma_W)
  
  # Hgamma here -> H in paper
  if (Zbin == 1) {
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  } else if (Zbin == 2) {
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  }

  # H here -> H_\delta in paper
  H = Hgamma[1:length(initd),1:length(initd)]
  
  # HI here -> Inverse of H_\delta in paper
  HI = ginv(H)
  
  # All derivatives with respect to a parameter of \delta AND \gamma.
  # Vargamma here -> H_\gamma in paper
  Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
  
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
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
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
      warning("This method used realV (for debugging)")
      newrow <- 2*realV[i]%*%XandW[i,]
      
      ##########################################################################
      mi = rbind(mi,newrow)
    }
    
    mi <- t(mi)
    
  } else if (Zbin == 2) {
    secder=t(-dlogis(XandW%*%gammaest))%*%prodvec
    
    # WM here -> M in paper
    WM = secder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-secder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    WM <- (2/n)*WM
    
    # WMI here -> M^{-1} in paper
    WMI = ginv(WM)
    
    diffvec = Z-plogis(XandW%*%gammaest)
    
    # mi here -> h_m(W, Z, \gamma^*) in paper
    mi = c()
    
    for(i in 1:n){
      newrow <- 2*diffvec[i,]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi <- t(mi)
  }
  
  # psii here -> \Psi in paper
  psii = -WMI%*%mi
  
  # gi here -> h_l(S_i, gamma, delta) in paper
  gi = c()
  
  for (i in 1:n) {
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
  
  #### Some useful variables, functions, etc. ####
  
  Y_trans_theta1 <- YJtrans(Y, parhat[totparl+4])
  Y_trans_theta2 <- YJtrans(Y, parhat[totparl+5])
  
  # Asymptotic linear representation of \delta
  Omega_delta <- - HI %*% (gi + Vargamma %*% psii)
  
  # Asymptotic linear representation of \gamma
  Omega_gamma <- psii
  
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
  
  ##############################################################################
  # For debugging: also create the 'true' FK, i.e. using the true parameters. Do
  # the same for GK and QK.
  
  Y_trans_theta1_true <- YJtrans(Y, parN[["sd"]][4])
  Y_trans_theta2_true <- YJtrans(Y, parN[["sd"]][5])
  
  # Compute the cumulative density of K according to our two-step model
  FK_true <- rep(0, n)
  for (j in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1_true[j] - M %*% parN[["beta"]])/parN[["sd"]][1]
    arg2_vector <- (Y_trans_theta2_true[j] - M %*% parN[["eta"]])/parN[["sd"]][2]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK_true[j] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                        pbivnorm(cbind(arg1_vector, arg2_vector), rho = parN[["sd"]][3]))/n
  }
  FK_true <- FK_true[order(FK_true)]
  
  par(mfrow = c(1, 1))
  plot(sort(Y), FK_true, type = 'l')
  lines(sort(Y), FK, type = "l", col = 'red')
  # Basically the same
  
  # Compute the censoring distribution of K. A ~ unif(-20, 20)
  GK_true <- pmax(0, sort(Y)/40 + 0.5)
  
  par(mfrow = c(1, 1))
  plot(surv_A)
  lines(sort(Y), 1 - GK_true, type = 'l', col = 'red')
  
  QK_true <- 1 - (1 - FK_true) * (1 - GK_true)
  
  plot(sort(Y), QK_true, type = 'l')
  lines(sort(Y), QK, type = 'l', col = 'red')
  # Some slightly larger differences but still not too much
  
  ##############################################################################

  #### Obtain expression for W_n(k) ####

  Omega_deltaT <- rbind(Omega_delta[c(totparl+4, 1:parl, totparl+1),], Omega_gamma)
  Omega_deltaC <- rbind(Omega_delta[c(totparl+5, (parl+1):(totparl), totparl+2),], Omega_gamma)
  Omega_deltarho <- Omega_delta[totparl+3,]
  
  
  tilde_Omega_delta <- rbind(Omega_deltaT, Omega_deltaC, Omega_deltarho)
  
  if (Zbin == 1) {
    partial.g.gamma <- -t(XandW)
  }

  Omega_vct <- function(k, dQ1k.bar, parl, totparl, parlgamma, n, partial.g.gamma) {
    
    Omegas <- rep(0, n)
    Omegas_KM <- rep(0, n)
    Omegas_2step <- rep(0, n)
    Omegas_EDF <- rep(0, n)
    
    Omegas_KM_true <- rep(0, n)
    
    
    ##### Precompute E[J(k; \gamma^*, \delta^*)] #####
    
    J <- matrix(nrow = (totparl + 5 + 2*parlgamma), ncol = n)
    
    for (i in 1:n) {
      tauT <- parhat[1] + as.matrix(X)[i,]%*%parhat[2:(parl-2)] + Z[i]*parhat[parl-1] +
        V[i]*parhat[parl]
      tauT <- as.numeric(tauT)
      
      tauC <- parhat[parl+1] + as.matrix(X)[i,]%*%parhat[(parl+2):(totparl-2)] +
        Z[i]*parhat[totparl-1] + V[i]*parhat[totparl]
      tauC <- as.numeric(tauC)
      
      JTi_tilde <- (1/parhat[totparl+1])*c(DthetaYJtrans(k, parhat[totparl + 4]),
                                           -XandW[i,1:(parl-2)],
                                           -Z[i],
                                           -V[i],
                                           -(Y_trans_theta1[i] - tauT)/parhat[totparl+1],
                                           -parhat[parl]*partial.g.gamma[,i])
      JTi <- dnorm((Y_trans_theta1[i] - tauT)/parhat[totparl+1])*JTi_tilde
      
      JCi_tilde <- (1/parhat[totparl+2])*c(DthetaYJtrans(k, parhat[totparl + 5]),
                                           -XandW[i,1:(parl-2)],
                                           -Z[i],
                                           -V[i],
                                           -(Y_trans_theta2[i] - tauC)/parhat[totparl+2],
                                           -parhat[totparl]*partial.g.gamma[,i])
      JCi <- dnorm((Y_trans_theta2[i] - tauC)/parhat[totparl+2])*JCi_tilde
      
      JTCi_tilde <- c(JTi_tilde, JCi_tilde, 1)
      JTCi <- dnorm2d((Y_trans_theta1[i] - tauT)/parhat[totparl+1],
                      (Y_trans_theta2[i] - tauC)/parhat[totparl+2],
                      rho = parhat[totparl + 3]) * JTCi_tilde
      
      J[,i] <- c(JTi, JCi, 0) - JTCi
    }
    
    EJ <- rowMeans(J)
    
    ##### Also precompute E[\Phi(...) + \Phi(...) - \Phi(..., ...)] #####
    
    # In the process, we will store intermediate results so we can use them
    # later on.
    Phi_T_k <- rep(0, n)
    Phi_C_k <- rep(0, n)
    Phi_TC_k <- rep(0, n)
    
    for (i in 1:n) {
      tauT <- parhat[1] + as.matrix(X)[i,]%*%parhat[2:(parl-2)] + Z[i]*parhat[parl-1] +
        V[i]*parhat[parl]
      tauT <- as.numeric(tauT)
      
      tauC <- parhat[parl+1] + as.matrix(X)[i,]%*%parhat[(parl+2):(totparl-2)] +
        Z[i]*parhat[totparl-1] + V[i]*parhat[totparl]
      tauC <- as.numeric(tauC)
      
      k_trans_theta1 <- YJtrans(k, parhat[totparl + 4])
      k_trans_theta2 <- YJtrans(k, parhat[totparl + 5])
      
      Phi_T_k[i] <- pnorm( (k_trans_theta1 - tauT)/parhat[totparl + 1])
      Phi_C_k[i] <- pnorm( (k_trans_theta2 - tauC)/parhat[totparl + 2])
      Phi_TC_k[i] <- pnorm2d((k_trans_theta1 - tauT)/parhat[totparl + 1],
                             (k_trans_theta2 - tauC)/parhat[totparl + 2],
                             rho = parhat[totparl + 3])
    }
    
    EPhi_k <- (1/n)*(sum(Phi_T_k) + sum(Phi_C_k) - sum(Phi_TC_k))
    
    ##### For each observation S_i, compute \Omega(k; S_i) #####
    
    for (i in 1:n) {
      
      ###### Computations for Omega_2step ######

      Omega_2step <- t(EJ) %*% tilde_Omega_delta[,i]
      
      ###### Computations for Omega_KM ######
      
      FK.k <- stepfun(sort(Y), c(0, FK))(k)
      
      Omega_KM <- (1 - FK.k)*(gk(min(k, Y[i]), Y, QK, dQ1k.bar) + 
                                ifelse((Y[i] <= k) & (Delta[i] + Xi[i] == 1),
                                       1/(1 - QK[which(sort(Y) == Y[i])]),
                                       0))
      
      ###### Computations for Omega_EDF ######
      
      Omega_EDF <- (Phi_T_k[i] + Phi_C_k[i] - Phi_TC_k[i]) - EPhi_k
      
      ##### Compute Omega #####
      Omegas[i] <- Omega_2step - Omega_KM + Omega_EDF
      
      Omegas_KM[i] <- Omega_KM
      Omegas_2step[i] <- Omega_2step
      Omegas_EDF[i] <- Omega_EDF
    }
    
    list(Omegas, Omegas_2step, Omegas_KM, Omegas_EDF)
  }
  
  #### Estimate variance of W_n(k) for a range of k ####
  
  # Set value for U (= infimum of right end points of event and censoring
  # distributions of K)
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
                                  .combine = 'rbind') %dopar% {
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
      
      # Compute estimate of variance of W_n(k). The computation of the mean is not
      # necessary since we know it will be zero theoretically.
      all_Omegas <- Omega_vct(k, dQ1k.bar, parl, totparl, parlgamma, n, partial.g.gamma)
      Omegas <- all_Omegas[[1]]
      
      # Since we know the mean will be zero, we can estimate the variance slightly more
      # accurately
      var.Omegas <- (1/n)*sum(Omegas^2)
      
      # For debugging
      var.Omegas_2step <- (1/n)*sum(all_Omegas[[2]]^2)
      var.Omegas_KM <- (1/n)*sum(all_Omegas[[3]]^2)
      var.Omegas_EDF <- (1/n)*sum(all_Omegas[[4]]^2)
      
      c(mean(Omegas), var.Omegas, var.Omegas_2step, var.Omegas_KM, var.Omegas_EDF,
        sum(Omegas))
    }
    
    mean_k <- mean_var_estimates[,1]
    var_k <- mean_var_estimates[,2]
    
    var_k.2step <- mean_var_estimates[,3]
    var_k.KM <- mean_var_estimates[,4]
    var_k.EDF <- mean_var_estimates[,5]
    
    sum_Omegas_k <- mean_var_estimates[,6]
    
  } else {
    
    stop("This version is not (yet) up-to-date with current progression.")
    
    mean_k <- rep(0, K_RANGE_SIZE)
    var_k <- rep(0, K_RANGE_SIZE)
    
    var.Omegas_2step <- rep(0, K_RANGE_SIZE)
    var.Omegas_KM <- rep(0, K_RANGE_SIZE)
    
    for (i in 1:K_RANGE_SIZE) {
      if ((i %% (K_RANGE_SIZE/10) == 0) & verbose) {
        message("Estimating variance: ", i/K_RANGE_SIZE*100, "% completion")
      }
      k <- k_range[i]
      
      # Compute Q1(k) and dQ1
      Q1k <- rep(0, n)
      for (j in 1:n) {
        Q1k[j] <- sum(Y < sort(Y)[j] & Delta + Xi == 1)/n
      }
      Q1k.bar <- 1 - Q1k
      
      dQ1k.bar <- rep(0, n)
      dQ1k.bar[1] <- 1 - Q1k.bar[1]
      dQ1k.bar[2:n] = Q1k.bar[2:n] - Q1k.bar[1:(n-1)]
      
      # Compute estimate of variance of W_n(k). The computation of the mean is not
      # necessary since we know it will be zero theoretically.
      all_Omegas <- Omega_vct(k, dQ1k.bar, parl, totparl, n, partial.g.gamma)
      Omegas <- all_Omegas[[1]]
      
      # Since we know the mean will be zero, we can estimate the variance slightly more
      # accurately
      mean_k[i] <- mean(Omegas)
      var_k[i] <- (1/length(Omegas))*sum(Omegas^2)
      
      # For debugging
      var.Omegas_2step[i] <- (1/length(all_Omegas[[2]]))*sum(all_Omegas[[2]]^2)
      var.Omegas_KM[i] <- (1/length(all_Omegas[[3]]))*sum(all_Omegas[[3]]^2)
    }
  }
  
  ##############################################################################
  # For debugging
  
  par(mfrow = c(1, 2))
  plot(k_range, mean_k)
  plot(k_range, var_k)
  
  par(mfrow = c(2, 2))
  
  plot(k_range, var_k)
  plot(k_range, var_k.2step)
  plot(k_range, var_k.KM)
  plot(k_range, var_k.EDF)
  
  par(mfrow = c(1, 1))
  
  # Would it be okay to ignore covariances?
  plot(k_range, var_k, ylim = c(0, 1), type = 'l', lwd = 2)
  lines(k_range, var_k.2step + var_k.KM + var_k.EDF, type = 'l', col = 'red',
        lwd = 2)
  legend("topleft", legend = c("total variance", "sum of comp. var."),
         col = c("black", "red"), lwd = 2)
  # Yes, approximately
  
  # Test if F_K(k) - F_KM(k) = (1/n) sum(\Omega)
  K_delta <- as.numeric((Delta == 1) | (Xi == 1)) # Censoring indicator for K = min(T,C)
  surv_km <- survfit(Surv(Y, K_delta) ~ 1, type = "kaplan-meier")
  cdf_km <- 1 - get_surv_prob(surv_km, sort(Y))
  
  par(mfrow = c(1, 1))
  plot(sort(Y), FK, type = 'l')
  lines(sort(Y), cdf_km, type = 's', col = "red")
  legend("topleft", legend = c("Estimated model", "Kaplan-Meier"),
         col = c("black", "red"), lty = 1)
  
  par(mfrow = c(1, 1))
  plot(sort(Y), (FK - cdf_km)^2, type = 'l', ylim = c(0, 10^(-3)))
  lines(k_range, sum_Omegas_k^2 / n, type = 'l', col = 'red')
  
  plot(k_range, sum_Omegas_k^2 / n, type = 'l', col = 'red')

  ##############################################################################
  
  #### Estimate distribution of LTCM ####
  
  # Compute dFK
  dFK = rep(0,n)
  dFK[1] = FK[1]
  dFK[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Compute LTCM
  LTCM_ESTIMATION_SIZE <- ltcm_estimation_size
  
  Wk <- matrix(nrow = LTCM_ESTIMATION_SIZE, ncol = n)
  for (i in 1:n) {
    # Since we do not compute the variance in each of the observed time points,
    # we will interpolate them.
    
    if (sort(Y)[i] > max(k_range)) {
      
      Wk_var <- var_k[K_RANGE_SIZE]
      
    } else if (sort(Y)[i] < min(k_range)) {
      
      Wk_var <- var_k[1]

    } else if (sort(Y)[i] %in% k_range) {
      
      Wk_var <- var_k[which(k_range == sort(Y)[i])]
      
    } else {
      
      # Find between which values of k the i'th ordered time lies
      k_range_prev_idx <- max(which(k_range <= sort(Y)[i]))
      k_range_next_idx <- min(which(k_range >= sort(Y)[i]))
      
      # Interpolate between estimated variances
      points <- data.frame(y = c(var_k[k_range_prev_idx], var_k[k_range_next_idx]),
                           x = c(k_range[k_range_prev_idx], k_range[k_range_next_idx]))
      
      # Make a prediction for the standard deviation at sort(Y)[i]
      Wk_var <- approx(points$x, points$y, xout = sort(Y)[i])$y
    }
    
    Wk[,i] <- rnorm(n = LTCM_ESTIMATION_SIZE, mean = 0, sd = sqrt(Wk_var))
  }
  
  LTCM_estimates <- rep(0, LTCM_ESTIMATION_SIZE)
  for (i in 1:LTCM_ESTIMATION_SIZE) {
    LTCM_estimates[i] <- sum((Wk[i,])^2 * weight_fun(sort(Y), U) * dFK)
  }
  
  #### Compute test statistic ####
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  K_delta <- as.numeric((Delta == 1) | (Xi == 1)) # Censoring indicator for K = min(T,C)
  surv_km <- survfit(Surv(Y, K_delta) ~ 1, type = "kaplan-meier")
  cdf_km <- 1 - get_surv_prob(surv_km, sort(Y))
  
  if (display.plot) {
    par(mfrow = c(1, 2))
    plot(sort(Y), FK, type = 'l', main = "Comparison of estimated model with KM",
         xlab = "K = min(T,C)", ylab = "F(K)")
    lines(sort(Y), cdf_km, type = 's', col = "red")
    legend("topleft", legend = c("Estimated model", "Kaplan-Meier"),
           col = c("black", "red"), lty = 1)
  }
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((cdf_km - FK)^2 * weight_fun(sort(Y), U) * w)
  
  if (display.plot) {
    hist(LTCM_estimates, main = "Histogram of LT_CM")
    abline(v = TCM, col = 'red')
  }
  
  #### Do inference ####
  
  pvalue <- sum(LTCM_estimates > TCM)/LTCM_ESTIMATION_SIZE
  q90 <- quantile(LTCM_estimates, 0.90)
  q95 <- quantile(LTCM_estimates, 0.95)
  
  list(TCM, LTCM_estimates, pvalue, q90, q95)
}


GOF_noBootstrap_EstimateRunDuration <- function(data, nruns, Zbin, Wbin, k_range_size,
                                                ltcm_estimation_size) {
  start.time <- Sys.time()
  out <- GOF_test_noBootstrap(data, Zbin, Wbin, k_range_size, ltcm_estimation_size, 
                              display.plot = FALSE)
  diff <- difftime(Sys.time(), start.time, units = "secs")
  diff <- as.integer(round(diff))
  duration <- round((diff*nruns)/60)
  if (duration > 60) {
    message("The simulation will take approximately ", floor(duration / 60),
            " hours and ", duration %% 60, " minutes")
  } else {
    message("The simulation will take approximately ", duration, " minutes")
  } 
}


GOF_SimulationNoBootstrap_typeIerror <- function(parN, nruns, iseed, Zbin, Wbin,
                                                 display.plot) {
  
  nbr_rejections_90 <- 0
  nbr_rejections_95 <- 0
  differences <- rep(0, nruns)
  for (run in 1:nruns) {
    message("\n" , "Starting run ", run, "\n")
    run_seed <- iseed + run
    
    data <- dat.sim.reg(n, parN, run_seed, Zbin, Wbin)
    
    out <- GOF_test_noBootstrap(data, Zbin, Wbin, k_range_size, ltcm_estimation_size, 
                                display.plot)
    
    if (out[[1]] > out[[4]]) {
      nbr_rejections_90 <- nbr_rejections_90 + 1
    }
    
    if (out[[1]] > out[[5]]) {
      nbr_rejections_95 <- nbr_rejections_95 + 1
    }
    
    differences[run] <- out[[1]] - mean(out[[2]])
  }
  
  filename <- paste0("typeIerror_iseed", iseed, "_nruns", nruns, ".csv")
  filename2 <- paste0("differences_iseed", iseed, "_nruns", nruns, ".csv")
  folder <- paste0("TypeIerror/TypeIerrorAsymptotic_n", n)
  
  # If folder doesn't exist, create it.
  if (!file.exists(folder)) {
    dir.create(folder)
  }
  
  write.csv(c(nbr_rejections_90/nruns, nbr_rejections_95/nruns), file = paste0(folder, "/", filename))
  write.csv(differences, file = paste0(folder, "/", filename2))
  
  if (display.plot) {
    par(mfrow = c(1, 1))
    hist(differences, main = "E[T_{CM} - T^*_{CM}]")
  }
  
  data.frame("prop. reject90" = nbr_rejections_90/nruns, "prop. reject95" = nbr_rejections_95/nruns)
}

################################################################################
#                           Deprecated functions                               #
################################################################################

# These functions are deprecated and no longer updated.
GOF_test <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  warning("GOF test function not using parallel computing is deprecated!")
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
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
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
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
  
  # Compute F_K(k; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FK <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
  }
  
  # Order the values
  FK <- FK[order(FK)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  K_delta <- as.numeric((Delta == 1) | (Xi == 1))
  surv_km <- survfit(Surv(Y, K_delta) ~ 1)
  cdf_km <- 1 - surv_km[["surv"]]
  
  # It can happen that the if there are two values of Y which are very close to
  # each other, the survfit function will view them as the same value and return
  # a vector of survival probabilities that is of length n - 1 or less.
  
  # If it happens...
  while (length(cdf_km) < length(FK)) {
    
    # Find where
    Y.sorted <- sort(Y)
    matched_until <- 0
    for (i in 1:length(surv_km[["time"]])) {
      if (Y.sorted[i] == surv_km[["time"]][i]) {
        matched_until <- i
      } else {
        break
      }
    }
    
    # The vectors will differ on index = 'matched_until' + 1. To fix this, we
    # just need to insert cdf_km[matched_until] on the (matched_until + 1)'th
    # position.
    cdf_km <- c(cdf_km[1:matched_until], cdf_km[matched_until:length(cdf_km)])
    surv_km[["time"]] <- c(surv_km[["time"]][1:matched_until], surv_km[["time"]][matched_until:length(surv_km[["time"]])])
    
    # The value of cdf_km increases by one each iteration, meaning that this 
    # loop will end eventually.
  }
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((cdf_km - FK)^2 * w)
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  TCMb_vector <- rep(0, B)
  for (b in 1:B) {
    
    # Create bootstrap data sample
    bootstrapseed <- B*iseed + b
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star <- rep(0, n)
    for (i in 1:n) {
      A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
    }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star==T.star)
    Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
    data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FK.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
      
      # Evaluate them in standard normal CDF, and compute the average
      FK.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
    }
    
    FK.star <- FK.star[order(FK.star)]
    
    w.star = rep(0,n)
    w.star[1] = FK.star[1]
    w.star[2:n] = FK.star[2:n]-FK.star[1:(n-1)]
    
    K_delta.star <- as.numeric((Delta.star == 1) | (Xi.star == 1))
    surv_km.star <- survfit(Surv(Y.star, K_delta.star) ~ 1)
    cdf_km.star <- 1 - surv_km.star[["surv"]]
    
    # Same issue as before
    while (length(cdf_km.star) < length(FK.star)) {
      
      Y.star.sorted <- sort(Y.star)
      matched_until <- 0
      for (i in 1:length(surv_km.star[["time"]])) {
        if (Y.star.sorted[i] == surv_km.star[["time"]][i]) {
          matched_until <- i
        } else {
          break
        }
      }
      
      surv_km.star[["time"]] <- c(surv_km.star[["time"]][1:matched_until], surv_km.star[["time"]][matched_until:length(surv_km.star[["time"]])])
      cdf_km.star <- c(cdf_km.star[1:matched_until], cdf_km.star[matched_until:length(cdf_km.star)])
    }
    
    TCMb_vector[b] <- n*sum((cdf_km.star - FK.star)^2 * w.star)
  }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}

# These functions use the first approach considered to obtain the Cramer -
# von Mises statistic. This approach will model F_Y instead of F_K.
# --> It was shown that they do not work and therefore shouldn't be used. These
#     functions will also no longer be updated.
GOF_test_newapproach <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
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
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
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
  
  # Compute F_Y(y; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FY <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FY[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
    
    # Take into account the administrative censoring. Note that, in theory, the
    # value for Y[i] should match exactly to some value in surv_A[["time"]], but
    # probably due to rounding, this is not always the case. Hence the necessity
    # for slightly complicated code below. The maximum with 1 is taken to keep
    # the code from throwing an error should, for some reason (maybe again due
    # to rounding), there would not be a value in surv_A[["time"]] which is
    # smaller than or equal to Y[i].
    FY[i] <- FY[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]) + 
      cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]
  }
  
  # Order the values
  FY <- FY[order(FY)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FY[1]
  w[2:n] = FY[2:n]-FY[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  ord <- 1:n
  EFY <- (ord-1)/n
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((EFY - FY)^2 * w)
  
  TCMb_vector <- rep(0, B)
  for (b in 1:B) {
    
    # Create bootstrap data sample
    bootstrapseed <- B*iseed + b
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star <- rep(0, n)
    for (i in 1:n) {
      A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
    }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star==T.star)
    Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
    data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FY.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
      
      # Evaluate them in standard normal CDF, and compute the average
      FY.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
      
      # Take into account the administrative censoring
      FY.star[i] <- FY.star[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]) + 
        cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]
    }
    
    FY.star <- FY.star[order(FY.star)]
    
    w.star = rep(0,n)
    w.star[1] = FY.star[1]
    w.star[2:n] = FY.star[2:n]-FY.star[1:(n-1)]
    
    ord <- 1:n
    EFY.star <- (ord-1)/n
    
    TCMb_vector[b] <- n*sum((EFY.star - FY.star)^2 * w.star)
  }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}
GOF_test_newapproach_parallel <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
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
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
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
  
  # Compute F_Y(y; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FY <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FY[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
    
    # Take into account the administrative censoring. Note that, in theory, the
    # value for Y[i] should match exactly to some value in surv_A[["time"]], but
    # probably due to rounding, this is not always the case. Hence the necessity
    # for slightly complicated code below. The maximum with 1 is taken to keep
    # the code from throwing an error should, for some reason (maybe again due
    # to rounding), there would not be a value in surv_A[["time"]] which is
    # smaller than or equal to Y[i].
    FY[i] <- FY[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]) + 
      cdf_A[max(c(which(surv_A[["time"]] <= Y[i])), 1)]
  }
  
  # Order the values
  FY <- FY[order(FY)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FY[1]
  w[2:n] = FY[2:n]-FY[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  ord <- 1:n
  EFY <- (ord-1)/n
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((EFY - FY)^2 * w)
  
  package.vector <- c("MASS", "nloptr", "VGAM", "pbivnorm", "survival")
  export.vector <- c("parl", "totparl")
  TCMb_vector <- foreach(b = 1:B,
                         .packages = package.vector,
                         .export = export.vector,
                         .combine = c) %dopar% {
                           
                           source(paste0(parent_directory, "/Functions_ad.R"))
                           source("Goodness-of-fit-test_functions.R")
                           
                           # Create bootstrap data sample
                           bootstrapseed <- B*iseed + b
                           set.seed(bootstrapseed)
                           
                           beta <- parhat[1:parl]
                           eta <- parhat[(parl+1):totparl]
                           sd <- parhat[(totparl+1):(totparl+5)]
                           gamma <- gammaest
                           
                           mu <- c(0,0)
                           sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
                           err <- mvrnorm(n, mu=mu , Sigma=sigma)
                           err1 = err[,1]
                           err2 = err[,2]
                           
                           T.star = IYJtrans(M %*% beta + err1, sd[4])
                           C.star = IYJtrans(M %*% eta + err2, sd[5])
                           A.star <- rep(0, n)
                           for (i in 1:n) {
                             A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
                           }
                           
                           Y.star = pmin(T.star,C.star,A.star)
                           Delta.star = as.numeric(Y.star==T.star)
                           Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
                           data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
                           
                           initd1 <- parhat
                           parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                                                 ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
                           
                           Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
                           Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
                           
                           # For each observed time y, compute F_K(k; gamma, delta)
                           FY.star <- rep(0, n)
                           
                           for (i in 1:n) {
                             
                             # Vector of arguments of first term and second term
                             arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
                             arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
                             
                             # Evaluate them in standard normal CDF, and compute the average
                             FY.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                                                 pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
                             
                             # Take into account the administrative censoring
                             FY.star[i] <- FY.star[i]*(1-cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]) + 
                               cdf_A[max(c(which(surv_A[["time"]] <= Y.star[i]), 1))]
                           }
                           
                           FY.star <- FY.star[order(FY.star)]
                           
                           w.star = rep(0,n)
                           w.star[1] = FY.star[1]
                           w.star[2:n] = FY.star[2:n]-FY.star[1:(n-1)]
                           
                           ord <- 1:n
                           EFY.star <- (ord-1)/n
                           
                           n*sum((EFY.star - FY.star)^2 * w.star)
                         }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}
GOF_SimulationTypeIerror_newapproach <- function(parN, nruns, B, iseed, Zbin, Wbin, parallel) {
  nbr_reject90 <- 0
  nbr_reject95 <- 0
  TCM_reject <- c()
  TCMb_reject <- c()
  
  # Set up for parallel computing if necessary
  if (parallel) {
    clust <- makeCluster(10)
    registerDoParallel(clust)
  }
  
  for (run in 1:nruns) {
    message(paste0("Currently busy with run ", run, " out of ", nruns, "."))
    
    data <- dat.sim.reg(n, parN, iseed + run, Zbin, Wbin)
    
    if (parallel) {
      out <- GOF_test_newapproach_parallel(data, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
    } else {
      out <- GOF_test_newapproach(data, B, iseed + run, Zbin, Wbin, display.plot = FALSE)
    }
    
    if (out[[3]]) {
      nbr_reject90 <- nbr_reject90 + 1
      TCM_reject <- cbind(TCM_reject, out[[1]])
      TCMb_reject <- cbind(TCMb_reject, out[[2]])
    }
    if (out[[4]]) {
      nbr_reject95 <- nbr_reject95 + 1
    }
  }
  
  # stop cluster after parallel computing
  if (parallel) {
    stopCluster(clust)
  }
  
  list(reject90 = nbr_reject90/nruns, reject95 = nbr_reject95/nruns,
       TCM_rejected = TCM_reject, TCMb_rejected = TCMb_reject)
}


# This function is exactly like GOF_test, but instead of bootstrap sampling A
# from its estimated distribution, we sample from its true distribution.
# --> This function is no longer used and will no longer be updated.
GOF_test_oracle <- function(data, B, iseed, Zbin, Wbin, display.plot) {
  
  if ((Zbin != 1 && Zbin != 2) && (Wbin != 1 && Wbin !=2) ) {
    stop("Invalid input")
  }
  
  # Estimate the parameter vectors gamma and delta
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
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
  } else {
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  }
  
  M = cbind(data[,4:(2+parl)],V)
  
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
  
  # Compute F_K(k; gamma, delta)
  
  Y_trans_theta1 <- YJtrans(Y, parhat[12])
  Y_trans_theta2 <- YJtrans(Y, parhat[13])
  
  # For each observed time y, compute F_K(k; gamma, delta)
  FK <- rep(0, n)
  
  for (i in 1:n) {
    
    # Vector of arguments of first term and second term
    arg1_vector <- (Y_trans_theta1[i] - M %*% parhat[1:4])/parhat[9]
    arg2_vector <- (Y_trans_theta2[i] - M %*% parhat[5:8])/parhat[10]
    
    # Evaluate them in standard normal CDF, and compute the average
    FK[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                   pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat[11]))/n
  }
  
  # Order the values
  FK <- FK[order(FK)]
  
  # Useful when computing the CM-statistic later on. In a sense, it can be seen to
  # replace dF_K.
  w = rep(0,n)
  w[1] = FK[1]
  w[2:n] = FK[2:n]-FK[1:(n-1)]
  
  # Also compute the cumulative distribution function based on KM-estimator
  K_delta <- as.numeric((Delta == 1) | (Xi == 1))
  surv_km <- survfit(Surv(Y, K_delta) ~ 1)
  cdf_km <- 1 - surv_km[["surv"]]
  
  # It can happen that the if there are two values of Y which are very close to
  # each other, the survfit function will view them as the same value and return
  # a vector of survival probabilities that is of length n - 1 or less.
  
  # If it happens...
  while (length(cdf_km) < length(FK)) {
    
    # Find where
    Y.sorted <- sort(Y)
    matched_until <- 0
    for (i in 1:length(surv_km[["time"]])) {
      if (Y.sorted[i] == surv_km[["time"]][i]) {
        matched_until <- i
      }
    }
    
    # The vectors will differ on index = 'matched_until' + 1. To fix this, we
    # just need to insert cdf_km[matched_until] on the (matched_until + 1)'th
    # position.
    cdf_km <- c(cdf_km[1:matched_until], cdf_km[matched_until:length(cdf_km)])
    surv_km[["time"]] <- c(surv_km[["time"]][1:matched_until], surv_km[["time"]][matched_until:length(surv_km[["time"]])])
    
    # The value of cdf_km increases by one each iteration, meaning that this 
    # loop will end eventually.
  }
  
  # Compute the Cramer - von Mises statistic
  TCM <- n*sum((cdf_km - FK)^2 * w)
  
  # In the following, an estimator of the distribution function of the adminis-
  # trative censoring times A will be useful.
  A_cens_ind <- as.numeric((Delta == 0) & (Xi == 0))
  surv_A <- survfit(Surv(Y, A_cens_ind) ~ 1)
  cdf_A = 1 - surv_A[["surv"]]
  
  TCMb_vector <- rep(0, B)
  for (b in 1:B) {
    
    # Create bootstrap data sample
    bootstrapseed <- B*iseed + b
    set.seed(bootstrapseed)
    
    beta <- parhat[1:parl]
    eta <- parhat[(parl+1):totparl]
    sd <- parhat[(totparl+1):(totparl+5)]
    gamma <- gammaest
    
    mu <- c(0,0)
    sigma <- matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
    err <- mvrnorm(n, mu=mu , Sigma=sigma)
    err1 = err[,1]
    err2 = err[,2]
    
    T.star = IYJtrans(M %*% beta + err1, sd[4])
    C.star = IYJtrans(M %*% eta + err2, sd[5])
    A.star = runif(n, 0, 8)
    # A.star <- rep(0, n)
    # for (i in 1:n) {
    #   A.star[i] <- surv_A[["time"]][max(which(cdf_A <= runif(1, 0, 1)))]
    # }
    
    Y.star = pmin(T.star,C.star,A.star)
    Delta.star = as.numeric(Y.star==T.star)
    Xi.star = ifelse(Y.star==T.star, 0, as.numeric(Y.star == C.star))
    data.star = cbind(Y.star, Delta.star, Xi.star, M, realV)
    
    initd1 <- parhat
    parhat.star <- nloptr(x0=initd1,eval_f=LikF,Y=Y.star,Delta=Delta.star,Xi=Xi.star,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),
                          ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Y.star_trans_theta1 <- YJtrans(Y.star, parhat.star[12])
    Y.star_trans_theta2 <- YJtrans(Y.star, parhat.star[13])
    
    # For each observed time y, compute F_K(k; gamma, delta)
    FK.star <- rep(0, n)
    
    for (i in 1:n) {
      
      # Vector of arguments of first term and second term
      arg1_vector <- (Y.star_trans_theta1[i] - M %*% parhat.star[1:4])/parhat.star[9]
      arg2_vector <- (Y.star_trans_theta2[i] - M %*% parhat.star[5:8])/parhat.star[10]
      
      # Evaluate them in standard normal CDF, and compute the average
      FK.star[i] <- sum(pnorm(arg1_vector) + pnorm(arg2_vector) -
                          pbivnorm(cbind(arg1_vector, arg2_vector), rho = parhat.star[11]))/n
    }
    
    FK.star <- FK.star[order(FK.star)]
    
    w.star = rep(0,n)
    w.star[1] = FK.star[1]
    w.star[2:n] = FK.star[2:n]-FK.star[1:(n-1)]
    
    K_delta.star <- as.numeric((Delta.star == 1) | (Xi.star == 1))
    surv_km.star <- survfit(Surv(Y.star, K_delta.star) ~ 1)
    cdf_km.star <- 1 - surv_km.star[["surv"]]
    
    # Same issue as before
    while (length(cdf_km.star) < length(FK.star)) {
      
      Y.star.sorted <- sort(Y.star)
      matched_until <- 0
      for (i in 1:length(surv_km.star[["time"]])) {
        if (Y.star.sorted[i] == surv_km.star[["time"]][i]) {
          matched_until <- i
        }
      }
      
      surv_km.star[["time"]] <- c(surv_km.star[["time"]][1:matched_until], surv_km.star[["time"]][matched_until:length(surv_km.star[["time"]])])
      cdf_km.star <- c(cdf_km.star[1:matched_until], cdf_km.star[matched_until:length(cdf_km.star)])
    }
    
    TCMb_vector[b] <- n*sum((cdf_km.star - FK.star)^2 * w.star)
  }
  
  if (display.plot) {
    hist(TCMb_vector, main = "Histogram of bootstrap Cramer-Von Mises statistics",
         xlab = c(),
         xlim = c(0, max(max(TCMb_vector) + 0.1, TCM)))
    abline(v = TCM, col = "red")
  }
  
  significant1 <- (TCM > quantile(TCMb_vector, prob=0.9))
  significant2 <- (TCM > quantile(TCMb_vector, prob=0.95))
  
  list(TCM = TCM, TCMb_vector = TCMb_vector, signif90 = significant1,
       signif95 = significant2)
}






