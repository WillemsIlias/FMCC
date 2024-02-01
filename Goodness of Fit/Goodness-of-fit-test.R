
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
parN = list(beta=c(2.5,2.6,1.8,2),eta=c(1.8,0.9,0.5,-2.2),sd=c(1.1,1.4,0.75,1, 0.5),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = parl - 1
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho","theta_1","theta_2")

n <- 2000
iseed <- 123
Zbin <- 1
Wbin <- 1
B <- 250
display.plot = TRUE


# Generate data and give it representative column names
data <- dat.sim.reg(n, parN, iseed, Zbin, Wbin)
colnames(data) <- c("Y", "delta", "xi", "intercept", "x1", "Z", "W", "realV")

################################################################################
#                    GOF test based parametric bootstrap                       #
################################################################################

# Get an idea of computing time of the full simulation. This function might take
# a while.
nruns <- 500
GOF_EstimateRunDuration(data, B, nruns, Zbin, Wbin, parallel = TRUE)

# Run the simulation for assessing the Type 1 error of the goodness-of-fit test.
# Each of the 5 parts to be ran will do 100 simulations.
nbr_parts_to_evaluate <- 5
nruns_per_part <- nruns/nbr_parts_to_evaluate

for (part in 3:nbr_parts_to_evaluate) {
  message("Running part ", part, " out of ", nbr_parts_to_evaluate, ".")
  message("")
  
  typeIerror_results <- GOF_SimulationTypeIerror(parN, nruns_per_part, B, iseed + (part - 1)*100,
                                                 Zbin, Wbin, parallel = TRUE,
                                                 save.bootstrap.distributions = TRUE)
  reject90 <- typeIerror_results[["reject90"]]
  reject95 <- typeIerror_results[["reject95"]]
  
  filename <- paste0("typeIerror_iseed", iseed + (part - 1)*nruns_per_part, "_",
                     nruns_per_part, "runs.csv")
  folder <- paste0("TypeIerror/TypeIerror_B", B, "_n", n)
  
  # If folder doesn't exist, create it.
  if (!file.exists(folder)) {
    dir.create(folder)
  }
  
  # Write results to disk
  write.csv(c(reject90, reject95), file = paste0(folder, "/", filename))

}

# Plot the l'th rejected bootstrap sample together with the bootstrap T_{CM}.
l <- 3
hist(typeIerror_results[["TCMb_rejected"]][,l],
     main = "Histogram of bootstrap TCM leading to rejection",
     xlab = c(),
     xlim = c(0, max(max(typeIerror_results[["TCMb_rejected"]][l]), typeIerror_results[["TCM_rejected"]][l])))
abline(v = typeIerror_results[["TCM_rejected"]][l], col = "red")


################################################################################
#             GOF test based on the limiting distribution of T_CM              #
################################################################################

# IMPORTANT NOTE:
# All of the methodology of this section is still in development and DOES NOT
# WORK PROPERLY. Use the bootstrap version of the goodness-of-fit test instead.

# To further prevent possible confusion, the GOF test based on the limiting 
# distribution will now throw an error when it is being called. For that reason,
# the code in this section no longer works.

# 
# Type-I error of test without bootstrap
#

n <- 1000
iseed <- 5656
Zbin <- 1
Wbin <- 1
k_range_size <- 100
ltcm_estimation_size <- 5000
display.plot <- TRUE
nruns <- 10

# Estimate duration of simulation
parN = list(beta=c(2.5, 2.6, 1.8, 2),
            eta=c(1.8, 0.9, 0.5, -2.2),
            sd=c(1.1, 1, 0.75, 1, 0.5),
            gamma=c(-1, 0.6, 2.3))

data <- dat.sim.reg.debug(n, parN, iseed, Zbin, Wbin)
colnames(data) <- c("Y", "delta", "xi", "intercept", "x1", "Z", "W", "realV")

clust <- makeCluster(10)
registerDoParallel(clust)

GOF_test_noBootstrap(data, Zbin, Wbin, k_range_size, ltcm_estimation_size,
                     display.plot = TRUE)

GOF_noBootstrap_EstimateRunDuration(data, nruns, Zbin, Wbin, k_range_size, ltcm_estimation_size)

# Run simulation
source("Goodness-of-fit-test_functions.R")
GOF_SimulationNoBootstrap_typeIerror(parN, nruns, iseed, Zbin, Wbin, display.plot)

stopCluster(clust)

#
# 'Quick and dirty' comparison of both methods
#

iseed <- 124
k_range_size <- 100
ltcm_estimation_size <- 1000
display.plot <- TRUE

# Regenerate the data set and compute the quantiles based on asymptotic theory.
data <- dat.sim.reg(n, parN, iseed, Zbin, Wbin)
colnames(data) <- c("Y", "delta", "xi", "intercept", "x1", "Z", "W", "realV")
out <- GOF_test_noBootstrap(data, Zbin, Wbin, k_range_size, ltcm_estimation_size, 
                            display.plot)
q90.as <- out[[4]]
q95.as <- out[[5]]

# Load the data set of bootstrapped T_CM and compute bootstrap quantiles
filename <- paste0(iseed, ".csv")
folder <- paste0("TypeIerror/BootstrapTCM_B", B, "_n", n)
bootstrap_tcm <- read.csv(paste0(folder, "/", filename))[,2]
q90.bs <- quantile(bootstrap_tcm, 0.90)
q95.bs <- quantile(bootstrap_tcm, 0.95)

comparison <- data.frame("q90 asymptotic" = q90.as, "q90 bootstrap" = q90.bs,
                         "q95 asymptotic" = q95.as, "q95 bootstrap" = q95.bs)
rownames(comparison) <- c("")
round(comparison, 3)

# Plot distributions and TCM
TCM <- out[[1]]
LTCM_distribution <- out[[2]]

par(mfrow = c(1, 2))
hist(bootstrap_tcm)
lines(abline(v = TCM, col = 'red'))
hist(LTCM_distribution)
lines(abline(v = TCM, col = 'red'))


################################################################################
#      GOF test for misspecified models: misspecified error distribution       #
################################################################################



# Copied from 'Misspecification_main.R'... >>>..................................

beta <- c(2.5, 2.6, 1.8, 2)
eta <- c(1.8,0.9,0.5,-2.2)
gamma <- c(-1, 0.6, 2.3)
sd <- c(1.1,1.4,0.75,1, 0.5)
iseed <- 12012000

parl <- length(beta)
totparl <- 2*parl
parlgamma <- (parl-1)

namescoef =  c("$\\beta_{T,0}$","$\\beta_{T,1}$","$\\alpha_T$","$\\lambda_T$",
               "$\\beta_{C,0}$","$\\beta_{C,1}$","$\\alpha_C$","$\\lambda_C$",
               "$\\sigma_T$","$\\sigma_C$","$\\rho$","$\\theta_1$","$\\theta_2$")

# ...

# 
# Skew-normal errors
#

# Set some more variables
psi.delta.theta <- c(-0.9, 0.97, 0.97, 1, 0.5)
par.skew_normal <- list(beta, eta, psi.delta.theta, gamma)

# ...

#
# Multivariate t distribution
#

par.df <- c(sd,3)
par.t <- list(beta, eta, par.df, gamma)

# ...

par.heteroscedastic <- list(beta, eta, sd, gamma)

# ........................<<<<<< Last line copied from 'Misspecification_main.R'



iseed <- 768266
type <- "heteroscedastic"
nruns <- 500
B <- 250
n <- 2000

nbr_parts_to_evaluate <- 5
nruns_per_part <- nruns/nbr_parts_to_evaluate
for (part in 1:nbr_parts_to_evaluate) {
  message("Running part ", part, " out of 5.")
  message("")
  
  if (type == "t") {
    par.vector <- par.t 
  } else if (type == "skew"){
    par.vector <- par.skew_normal
  } else {
    par.vector <- par.heteroscedastic
  }
  misspec_results <- 
    GOF_SimulationMisspecification(type, par.vector, nruns_per_part, B, iseed + (part - 1)*nruns_per_part, 
                                   Zbin, Wbin, parallel = TRUE)
  reject90 <- misspec_results[["reject90"]]
  reject95 <- misspec_results[["reject95"]]
  
  filename <- paste0("misspec_", type, "_iseed", iseed + (part - 1)*nruns_per_part, "_",
                     nruns_per_part, "runs_B", B, ".csv")
  
  # put in proper folder
  if (type == "heteroscedastic") {
    folder <- paste0("Power/Misspec_hs_B", B, "_n", n)
  } else {
    folder <- paste0("Power/Misspec_", type, "_B", B, "_n", n)
  }
  
  # If folder doesn't exist, create it.
  if (!file.exists(folder)) {
    dir.create(folder)
  }
  
  write.csv(c(reject90, reject95), file = paste0(folder, "/", filename))
  
}


################################################################################
#       GOF test for misspecified models: misspecified control function        #
################################################################################

source("Goodness-of-fit-test_functions.R")

beta <- c(2.5, 2.6, 1.8, 2)
eta <- c(1.8,0.9,0.5,-2.2)
gamma <- c(-1, 0.6, 2.3)
sd <- c(1.1,1.4,0.75,1, 0.5)

parl <- length(beta)
totparl <- 2*parl
parlgamma <- (parl-1)

namescoef =  c("$\\beta_{T,0}$","$\\beta_{T,1}$","$\\alpha_T$","$\\lambda_T$",
               "$\\beta_{C,0}$","$\\beta_{C,1}$","$\\alpha_C$","$\\lambda_C$",
               "$\\sigma_T$","$\\sigma_C$","$\\rho$","$\\theta_1$","$\\theta_2$")

iseed <- 768266
control_function <- "probit"
nruns <- 500
Zbin <- 2
Wbin <- 1
B <- 250
n <- 2000

nbr_parts_to_evaluate <- 5
nruns_per_part <- nruns/nbr_parts_to_evaluate
for (part in 4:nbr_parts_to_evaluate) {
  message("Running part ", part, " out of 5.")
  message("")
  
  if ((control_function == "probit") | (control_function == "cloglog")) {
    par.vector <- list(beta, eta, sd, gamma)
  } else {
    stop("invalid control function specified")
  }
  misspec_results <- 
    GOF_SimulationMisspecification(control_function, par.vector, nruns_per_part, B, iseed + (part - 1)*nruns_per_part, 
                                   Zbin, Wbin, parallel = TRUE)
  reject90 <- misspec_results[["reject90"]]
  reject95 <- misspec_results[["reject95"]]
  
  filename <- paste0("misspec_", control_function, "_iseed", iseed + (part - 1)*nruns_per_part, "_",
                     nruns_per_part, "runs_B", B, ".csv")
  
  # put in proper folder
  folder <- paste0("Power/Misspec_", control_function, "_B", B, "_n", n)
  
  # If folder doesn't exist, create it.
  if (!file.exists(folder)) {
    dir.create(folder)
  }
  
  write.csv(c(reject90, reject95), file = paste0(folder, "/", filename))
  
}











