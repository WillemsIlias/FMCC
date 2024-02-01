# the original data can be found at: https://www.upjohn.org/data-tools/employment-research-data-center/national-jtpa-study
# The function used to clean the data can be found at: https://github.com/GillesCrommen/DCC

# Clear workspace
rm(list = ls())

# Load necessary packages, functions and data
fulldata<-read.csv("Data/clean_dataset_JTPA_ad.csv")
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

# Select all white, married men who do not have children
dataset = fulldata[(fulldata$children == 0 &
                    fulldata$married  == 1 &
                    fulldata$male     == 1 &
                    fulldata$white    == 1), -c(1,4,5,7,8)]

# Create vectors for the observed times and censoring indicators
Y <- as.matrix(log(dataset$days))
Delta <- as.matrix(dataset$delta)
Xi <- as.matrix(dataset$xi)

# Create a vector that represents the intercept
dataset$intercept <- rep(1,nrow(dataset))

# Create data matrix X of unconfounded predictors. The columns of X are:
# [Intercept, age, hsged]
X <- as.matrix(subset(dataset, select = c(ncol(dataset),1,2)))

# Define some useful variables
parl <- ncol(X)+2
totparl <- 2*parl
parlgamma <- parl-1
n <- nrow(dataset)

# Create data matrix of confounded variable. Z is the indicator whether the
# participant actually participated in the study. (0 for no participation,
# 1 for participation). 
Z <- as.matrix(dataset$jtpa)

# Create data matrix of instrumental variable. W is the indicator whether the 
# participant was in the control (W = 0) or treatment group (W = 1).
W <- as.matrix(dataset$treatment)

# Create data matrix
XandW <- as.matrix(cbind(X,W))
data <- as.matrix(cbind(Y, Delta, Xi, X, Z, W))
namescoef <- c("beta_{T,0}","beta_{T,1}","beta_{T,2}",
               "alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","beta_{C,2}",
               "alpha_C","lambda_C","sigma_T",
               "sigma_C","rho","theta_1","theta_2")


# Fitting logistic regression model
fit_log <- glm(Delta ~ X[,2] + X[,3], family = "binomial")
summary(fit_log)

# We follow the approach as in Chapter 5 of Crommen, Van Keilegom (2022). Hence,
# we will be comparing the models (1) assuming independence (2) assuming
# dependence but no confounding and (3) our model. It would not make sense to
# compare our model with the model of Crommen, Van Keilegom (2022) directly.
# However, we will compare the estimated survival curves.

init.value.theta_1 <- 1
init.value.theta_2 <- 2
DataApplicationJPTA(data, init.value.theta_1, init.value.theta_2) # Takes about 1 minute to run

# Goodness-of-fit test

# Load GOF test functions
source("Goodness of Fit/Goodness-of-fit-test_functions.R")

# Reorder the colums in the dataset
data_GOF <- dataset
data_GOF <- as.matrix(data_GOF[,c(7, 5, 6, 8, 1:4)])

# Define the initial seed
iseed <- 35497438

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
GOF_test_parallel(data_GOF, B, iseed, Zbin, Wbin, display.plot = TRUE)
stopCluster(clust)




