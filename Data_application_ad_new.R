# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
# fulldata <- read.csv("hiip_data_041515.csv")
library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
library(VGAM)
library(doParallel)
library(foreach)
library(pbivnorm)
library(survival)
source("Functions_ad.r")
source("Goodness of Fit/Goodness-of-fit-test_functions.R")

# Preprocess the data
dataset <- fulldata[(fulldata$followup_days != "M" & fulldata$followup_days != "I"),]
dataset$followup_days <- as.numeric(dataset$followup_days)
dataset <- dataset[dataset$followup_days > 0, ]
dataset <- dataset[1:50000,]
dataset$age_entry <- (dataset$age_entry - mean(dataset$age_entry))/sd(dataset$age_entry)
Y <- as.matrix(log(dataset$followup_days))
Delta <- as.matrix(dataset$final_vital_stat)
d1 <- as.numeric(dataset$final_vital_stat == 1) #death from breast cancer
d2 <- as.numeric (dataset$final_vital_stat == 2) #death from other cause
da <- as.numeric (dataset$final_vital_stat == 0) # alive or lost to follow up
dataset$intercept <- rep(1,nrow(dataset))

Z <- as.matrix(dataset$cohort == 1) # selected + participated
W <- as.matrix(dataset$cohort != 3) # selected
X <- as.matrix(subset(dataset, select = c(ncol(dataset), age_entry)))

XandW <- as.matrix(cbind(X,W))
n <- nrow(dataset)
data <- as.matrix(cbind(Y,d1,d2,X,Z,W))

namescoef <- c("beta_{T1,0}","beta_{T1,1}", "alpha_{T1}","lambda_{T1}",
               "beta_{T2,0}","beta_{T2,1}", "alpha_{T2}","lambda_{T2}",
               "sigma_{T1}","sigma_{T2}","rho","theta_1","theta_2")

# Define some useful variables
parl <- ncol(X) + 2
totparl <- 2*parl
parlgamma <- parl - 2
n <- nrow(data)

# Define initial values for the transformation parameters
init.value.theta_1 <- 1
init.value.theta_2 <- 1

# For testing
if (FALSE) {
  data <- read.csv("test_data_set.csv")
  data <- data[1:300,]
  parl <- 4
  totparl <- 8
  parlgamma <- 3
  n <- nrow(data)
}

# Run the data application
DataApplicationBC(data, init.value.theta_1, init.value.theta_2,
                  multiple.starting.points = TRUE)

# Goodness-of-fit test

# Put data in matrix format
data.matrix <- as.matrix(data)

# Define the initial seed
iseed <- 387439821

# Variables indicating whether Z and W are binary or continuous variables. In
# this case, both are binary
Zbin <- 2
Wbin <- 2

# Number of bootstrap samples to be used in the bootstrap procedure for the
# goodness-of-fit test
B <- 500

# Necessary for proper execution of "GOF_test_parallel(...)"
parent_directory <- getwd()

# Run the goodness-of-fit test
clust <- makeCluster(10)
registerDoParallel(clust)
GOF_test_parallel(data.matrix, B, iseed, Zbin, Wbin, display.plot = TRUE)
stopCluster(clust)

