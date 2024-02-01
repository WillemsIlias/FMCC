
################# Yeo-Johnson transformation functions #########################

log_transform = function(y)   
{ 
  transform_y = (y>0)*y+(y<=0)*1 # sy>0 --> condition is a check
  return(log(transform_y)) 
} 

power_transform = function(y,pw) 
{ 
  transform_y = (y>0)*y+(y<=0)*1  # sy>0
  return(transform_y^pw) 
} 

YJtrans = function(y,theta) # Yeo-Johnson transformation 
{ 
  sg = y>=0 
  if (theta==0) {temp = log_transform(y+1)*sg+(1-sg)*(0.5-0.5*(y-1)^2)} 
  if (theta==2) {temp = sg*(-0.5+0.5*(y+1)^2)-log_transform(-y+1)*(1-sg)} 
  if ((theta!=0) & (theta!=2)) {temp = sg*(power_transform(y+1,theta)-1)/theta+(1-sg)*(1-power_transform(-y+1,2-theta))/(2-theta)} 
  return(temp) 
} 

IYJtrans = function(y,theta) # Inverse of Yeo-Johnson transformation 
{ 
  sg = y>=0 
  if (theta==0) {temp =(exp(y)-1)*sg+(1-sg)*(1-power_transform(-2*y+1,0.5))} 
  if (theta==2) {temp = sg*(-1+power_transform(2*y+1,0.5))+(1-exp(-y))*(1-sg)} 
  if ((theta!=0) & (theta!=2)) {temp = sg*(power_transform(abs(theta)*y+1,1/theta)-1)+(1-sg)*(1-power_transform(1-(2-theta)*y,1/(2-theta)))} 
  return(temp) 
} 

DYJtrans = function(y,theta) # Derivative of Yeo-Johnson transformation 
{ 
  sg = y>=0 
  temp = power_transform(y+1,theta-1)*sg+power_transform(-y+1,1-theta)*(1-sg) 
  return(temp) 
} 


######################## Integrals CIF calculation #############################

# competing risk 1
integral1 <- function(x, theta1, theta2, sigma1, sigma2, rho12,beta1,beta2,M){
  transY.T1=YJtrans(x,theta1)
  DtransY.T1=DYJtrans(x,theta1)
  transY.T2=YJtrans(x,theta2)
  DtransY.T2=DYJtrans(x,theta2)
  
  b_T1 = (transY.T1-(M%*%beta1))/sigma1 # b_T1
  b_T2 = (transY.T2-(M%*%beta2))/sigma2
  
  f_T1 = 1/sigma1*dnorm(b_T1)*DtransY.T1*pnorm(q=-(transY.T2-M%*%beta2-rho12*sigma2*b_T1)/(sigma2*(1-rho12^2)^(1/2)))
  return(f_T1)
}

# competing risk 2
integral2 <- function(x, theta1, theta2, sigma1, sigma2, rho12,beta1,beta2,M){
  transY.T1=YJtrans(x,theta1)
  DtransY.T1=DYJtrans(x,theta1)
  transY.T2=YJtrans(x,theta2)
  DtransY.T2=DYJtrans(x,theta2)
  
  
  b_T1 = (transY.T1-(M%*%beta1))/sigma1 # b_T1
  b_T2 = (transY.T2-(M%*%beta2))/sigma2
  
  f_T2 = 1/sigma2*dnorm(b_T2)*DtransY.T2*pnorm(q=-(transY.T1-M%*%beta1-rho12*sigma1*b_T2)/(sigma1*(1-rho12^2)^(1/2)))
  return(f_T2)
}


###################### data simulating function ################################

dat.sim.reg = function(n,par,iseed,Zbin,Wbin){
  
  
  set.seed(iseed)
  beta1 = par[[1]]
  beta2 = par[[2]]
  beta3 = par[[3]]
  sd = par[[4]]
  gamma = par[[5]]
  
  # bivariate normal distribution of error terms
  mu = c(0,0,0)
  sigma = matrix(c(sd[1]^2,sd[1]*sd[2]*sd[4],sd[1]*sd[3]*sd[5], sd[1]*sd[2]*sd[4],sd[2]^2,sd[2]*sd[3]*sd[6],sd[1]*sd[3]*sd[5],sd[2]*sd[3]*sd[6], sd[3]^2),ncol=3)
  err = mvrnorm(n, mu =mu , Sigma=sigma)
  
  # error T and error C
  err1 = err[,1]
  err2 = err[,2]
  err3 = err[,3]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = sample(c(0,1), n, replace = TRUE)
  
  if (Wbin==2) { # Bernoulli with p =0.5
    W = sample(c(0,1), n, replace = TRUE) # sample 0 and 1 with equal probability
  } else if (Wbin==1) {
    W = runif(n,0,2) #Uniform[0,2]
  }
  
  XandW=as.matrix(cbind(x0,x1,W)) # W vector
  
  if (Zbin==2) {  # Z is binary
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  } else if (Zbin==1) {# Z is continuous
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n)  # matrix containing all covariates
  T1 = IYJtrans(Mgen%*%beta1+err1,sd[7]) 
  T2 = IYJtrans(Mgen%*%beta2+err2,sd[8]) 
  T3 = IYJtrans(Mgen%*%beta3+err3,sd[9])
  A = runif(n,0,15) # administrative censoring
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix

  Y = pmin(T1,T2,T3,A) # observed non-transformed time
  d1 = as.numeric(Y==T1) # censoring indicator
  d2 = ifelse(d1==0,as.numeric(Y==T2),0)
  d3 = ifelse(d1==1 | d2==1,0,as.numeric(Y==T3))
  da = 1-d1-d2-d3
  data = cbind(Y,d1,d2,d3,da,M,realV) # data consisting of observed time,
  # censoring indicators, all data and the control function
  
  return(data)
}


######################## likelihood function ###################################


# maximum likelihood for Gamma (Z continuous)
LikGamma1 = function(par,Y,M){ 
  W=as.matrix(M)
  gamma= as.matrix(par)
  
  tot = (Y-W%*%gamma)^2
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(Logn)
}

# maximum likelihood for Gamma (Z discrete)
LikGamma2 = function(par,Y,M){ 
  W=as.matrix(M)
  gamma= as.matrix(par)
  
  tot = (plogis(W%*%gamma)^Y)*((1-plogis(W%*%gamma))^(1-Y))
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}


###########################################################
# joint model with dependent censoring and transformation #
###########################################################

# - Used in second step of 2-step estimation
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

LikF.cmprsk = function(par,Y,d1,d2,d3,da,M){ 
  M=as.matrix(M)
  k = ncol(M)
  l = 2*k
  p = 3*k
  v = k+1
  w = 2*k+1
  beta1 = as.matrix(par[1:k])
  beta2 = as.matrix(par[v:l])
  beta3 = as.matrix(par[w:p])
  sigma1 = par[p+1]
  sigma2 = par[p+2]
  sigma3 = par[p+3]
  rho12 = par[p+4]
  rho13 = par[p+5]
  rho23 = par[p+6]
  sigma <- matrix(c(sigma1^2,sigma1*sigma2*rho12,sigma1*sigma3*rho13, sigma1*sigma2*rho12,sigma2^2,sigma2*sigma3*rho23,sigma1*sigma3*rho13,sigma2*sigma3*rho23,sigma3^2),ncol=3)
  theta_1 = par[p+7]
  theta_2 = par[p+8]
  theta_3 = par[p+9]
  
  if (is.positive.definite(sigma,tol=1e-30)){ # Not every sigma, rho combination results in valid variance covariance matrix (check positive definite)
    transY.T1=YJtrans(Y,theta_1)
    DtransY.T1=DYJtrans(Y,theta_1)
    transY.T2=YJtrans(Y,theta_2)
    DtransY.T2=DYJtrans(Y,theta_2)
    transY.T3=YJtrans(Y,theta_3)
    DtransY.T3=DYJtrans(Y,theta_3)
    
    b_T1 = (transY.T1-(M%*%beta1))/sigma1 # b_T1
    b_T2 = (transY.T2-(M%*%beta2))/sigma2 
    b_T3 = (transY.T3-(M%*%beta3))/sigma3
    
    rho1 = (rho23-rho12*rho13)/((1-rho12^2)*(1-rho13^2))^(1/2)
    rho2 = (rho13-rho12*rho23)/((1-rho12^2)*(1-rho23^2))^(1/2)
    rho3 = (rho12-rho23*rho13)/((1-rho23^2)*(1-rho13^2))^(1/2)
    
    f_T1 = 1/sigma1*dnorm(b_T1)*DtransY.T1*pbinorm(q1=-(transY.T2-M%*%beta2-rho12*sigma2*b_T1)/(sigma2*(1-rho12^2)^(1/2)),q2=-(transY.T3-M%*%beta3-rho13*sigma3*b_T1)/(sigma3*(1-rho13^2)^(1/2)),cov12=rho1)
    f_T2 = 1/sigma2*dnorm(b_T2)*DtransY.T2*pbinorm(q1=-(transY.T1-M%*%beta1-rho12*sigma1*b_T2)/(sigma1*(1-rho12^2)^(1/2)),q2=-(transY.T3-M%*%beta3-rho23*sigma3*b_T2)/(sigma3*(1-rho23^2)^(1/2)),cov12=rho2)
    f_T3 = 1/sigma3*dnorm(b_T3)*DtransY.T3*pbinorm(q1=-(transY.T1-M%*%beta1-rho13*sigma1*b_T3)/(sigma1*(1-rho13^2)^(1/2)),q2=-(transY.T2-M%*%beta2-rho23*sigma2*b_T3)/(sigma2*(1-rho23^2)^(1/2)),cov12=rho3)
    
    f_Ta = pmnorm(cbind(-b_T1, -b_T2, -b_T3),mean=c(0,0,0), varcov=sigma )
    
    tot = f_T1^d1*f_T2^d2*f_T3^d3*f_Ta^da
    p1 = pmax(tot,1e-100)   
    Logn = sum(log(p1)); 
    return(-Logn)
  }else{
    return(2^{31}-1) # if not positive definite --> return a very large value
  }
  
}
################################
# Independent model assumption #
################################

# - Does assume endogeneity
# - Doesn't assume dependent censoring or dependency between competing risks
# - Uses Yeo-Johnson transformation

LikI.bis = function(par,Y,d1,d2,d3,da,M){ 
  M=as.matrix(M)
  k = ncol(M)
  l = 2*k
  p = 3*k
  v = k+1
  w = 2*k+1
  beta1 = as.matrix(par[1:k])
  beta2 = as.matrix(par[v:l])
  beta3 = as.matrix(par[w:p])
  sigma1 = par[p+1]
  sigma2 = par[p+2]
  sigma3 = par[p+3]
  rho12 = 0 # all correlations are assumed equal to zero
  rho13 = 0
  rho23 = 0
  sigma <- matrix(c(sigma1^2,sigma1*sigma2*rho12,sigma1*sigma3*rho13, sigma1*sigma2*rho12,sigma2^2,sigma2*sigma3*rho23,sigma1*sigma3*rho13,sigma2*sigma3*rho23,sigma3^2),ncol=3)
  theta_1 = par[p+4]
  theta_2 = par[p+5]
  theta_3 = par[p+6]
  
  transY.T1=YJtrans(Y,theta_1)
  DtransY.T1=DYJtrans(Y,theta_1)
  transY.T2=YJtrans(Y,theta_2)
  DtransY.T2=DYJtrans(Y,theta_2)
  transY.T3=YJtrans(Y,theta_3)
  DtransY.T3=DYJtrans(Y,theta_3)
  
  b_T1 = (transY.T1-(M%*%beta1))/sigma1 # b_T1
  b_T2 = (transY.T2-(M%*%beta2))/sigma2 
  b_T3 = (transY.T3-(M%*%beta3))/sigma3
  
  rho1 = (rho23-rho12*rho13)/((1-rho12^2)*(1-rho13^2))^(1/2)
  rho2 = (rho13-rho12*rho23)/((1-rho12^2)*(1-rho23^2))^(1/2)
  rho3 = (rho12-rho23*rho13)/((1-rho23^2)*(1-rho13^2))^(1/2)
  
  f_T1 = 1/sigma1*dnorm(b_T1)*DtransY.T1*pbinorm(q1=-(transY.T2-M%*%beta2-rho12*sigma2*b_T1)/(sigma2*(1-rho12^2)^(1/2)),q2=-(transY.T3-M%*%beta3-rho13*sigma3*b_T1)/(sigma3*(1-rho13^2)^(1/2)),cov12=rho1)
  f_T2 = 1/sigma2*dnorm(b_T2)*DtransY.T2*pbinorm(q1=-(transY.T1-M%*%beta1-rho12*sigma1*b_T2)/(sigma1*(1-rho12^2)^(1/2)),q2=-(transY.T3-M%*%beta3-rho23*sigma3*b_T2)/(sigma3*(1-rho23^2)^(1/2)),cov12=rho2)
  f_T3 = 1/sigma3*dnorm(b_T3)*DtransY.T3*pbinorm(q1=-(transY.T1-M%*%beta1-rho13*sigma1*b_T3)/(sigma1*(1-rho13^2)^(1/2)),q2=-(transY.T2-M%*%beta2-rho23*sigma2*b_T3)/(sigma2*(1-rho23^2)^(1/2)),cov12=rho3)
  
  f_Ta = pmnorm(cbind(-b_T1, -b_T2, -b_T3),mean=c(0,0,0), varcov=sigma )
  
  tot = f_T1^d1*f_T2^d2*f_T3^d3*f_Ta^da
  p1 = pmax(tot,1e-100) 
  Logn = sum(log(p1)); # calculate loglikelihood
  return(-Logn)
  
  
}

# - Does assume endogeneity
# - Doesn't assume dependent censoring but does assume dependency between competing risks
# - Uses Yeo-Johnson transformation

LikI.cmprsk = function(par,Y,d1,d2,d3,da,M){ 
  M=as.matrix(M)
  k = ncol(M)
  l = 2*k
  p = 3*k
  v = k+1
  w = 2*k+1
  beta1 = as.matrix(par[1:k])
  beta2 = as.matrix(par[v:l])
  beta3 = as.matrix(par[w:p])
  sigma1 = par[p+1]
  sigma2 = par[p+2]
  sigma3 = par[p+3]
  rho12 = par[p+4] # correlation between competing risks is non-zero
  rho13 = 0 # independent censoring is assumed
  rho23 = 0
  sigma <- matrix(c(sigma1^2,sigma1*sigma2*rho12,sigma1*sigma3*rho13, sigma1*sigma2*rho12,sigma2^2,sigma2*sigma3*rho23,sigma1*sigma3*rho13,sigma2*sigma3*rho23,sigma3^2),ncol=3)
  theta_1 = par[p+5]
  theta_2 = par[p+6]
  theta_3 = par[p+7]
  
  if (is.positive.definite(sigma,tol=1e-30)){
    
    transY.T1=YJtrans(Y,theta_1)
    DtransY.T1=DYJtrans(Y,theta_1)
    transY.T2=YJtrans(Y,theta_2)
    DtransY.T2=DYJtrans(Y,theta_2)
    transY.T3=YJtrans(Y,theta_3)
    DtransY.T3=DYJtrans(Y,theta_3)
    
    b_T1 = (transY.T1-(M%*%beta1))/sigma1 # b_T1
    b_T2 = (transY.T2-(M%*%beta2))/sigma2 
    b_T3 = (transY.T3-(M%*%beta3))/sigma3
    
    rho1 = (rho23-rho12*rho13)/((1-rho12^2)*(1-rho13^2))^(1/2)
    rho2 = (rho13-rho12*rho23)/((1-rho12^2)*(1-rho23^2))^(1/2)
    rho3 = (rho12-rho23*rho13)/((1-rho23^2)*(1-rho13^2))^(1/2)
    
    f_T1 = 1/sigma1*dnorm(b_T1)*DtransY.T1*pbinorm(q1=-(transY.T2-M%*%beta2-rho12*sigma2*b_T1)/(sigma2*(1-rho12^2)^(1/2)),q2=-(transY.T3-M%*%beta3-rho13*sigma3*b_T1)/(sigma3*(1-rho13^2)^(1/2)),cov12=rho1)
    f_T2 = 1/sigma2*dnorm(b_T2)*DtransY.T2*pbinorm(q1=-(transY.T1-M%*%beta1-rho12*sigma1*b_T2)/(sigma1*(1-rho12^2)^(1/2)),q2=-(transY.T3-M%*%beta3-rho23*sigma3*b_T2)/(sigma3*(1-rho23^2)^(1/2)),cov12=rho2)
    f_T3 = 1/sigma3*dnorm(b_T3)*DtransY.T3*pbinorm(q1=-(transY.T1-M%*%beta1-rho13*sigma1*b_T3)/(sigma1*(1-rho13^2)^(1/2)),q2=-(transY.T2-M%*%beta2-rho23*sigma2*b_T3)/(sigma2*(1-rho23^2)^(1/2)),cov12=rho3)
    
    f_Ta = pmnorm(cbind(-b_T1, -b_T2, -b_T3),mean=c(0,0,0), varcov=sigma )
    
    tot = f_T1^d1*f_T2^d2*f_T3^d3*f_Ta^da
    p1 = pmax(tot,1e-100) 
    Logn = sum(log(p1)); # calculate loglikelihood
    return(-Logn)
  }else{
    return(2^{31}-1)
  }
}

####################
# Data application #
####################

# Likelihood functions data application (V=0 method)

LikFG.app = function(par,Y,Delta,Xi,M){ 
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  rho = par[l+5]
  theta_1 = par[l+6]
  theta_2 = par[l+7]
  gamma = as.matrix(par[(l+8):(l+7+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=M[,k+2]
  Vest <- rep(0,length(Z))
  Vest[W==1] <-(1-Z[W==1])*((1+exp(X[W==1,]%*%gamma))*log(1+exp(X[W==1,]%*%gamma))-(X[W==1,]%*%gamma)*exp(X[W==1,]%*%gamma))-Z[W==1]*((1+exp(-(X[W==1,]%*%gamma)))*log(1+exp(-(X[W==1,]%*%gamma)))+(X[W==1,]%*%gamma)*exp(-(X[W==1,]%*%gamma)))
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((transY.C-rho*sigma2/sigma1*transY.T)-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (transY.C-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((transY.T-rho*sigma1/sigma2*transY.C)-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*(((1/sigma2)*dnorm(z3)*(1-pnorm(z4))*DtransY.C)^Xi)*((pbinorm(q1=-z1,q2=-z3,cov12=rho))^(1-(Delta+Xi))) # likelihood
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}


LikIGamma.app = function(par,Y,Delta,Xi,M){ 
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  theta_1 = par[l+5]
  theta_2 = par[l+6]
  gamma = as.matrix(par[(l+7):(l+6+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=M[,k+2]
  Vest <- rep(0,length(Z))
  Vest[W==1] <-(1-Z[W==1])*((1+exp(X[W==1,]%*%gamma))*log(1+exp(X[W==1,]%*%gamma))-(X[W==1,]%*%gamma)*exp(X[W==1,]%*%gamma))-Z[W==1]*((1+exp(-(X[W==1,]%*%gamma)))*log(1+exp(-(X[W==1,]%*%gamma)))+(X[W==1,]%*%gamma)*exp(-(X[W==1,]%*%gamma)))
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = (transY.C-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1))*DtransY.C)^Xi)*((pbinorm(q1=-z1,q2=-z2,cov12=0))^(1-(Delta+Xi)))
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}



######################## simulation function ###################################
SimulationCI22_Sara = function(n, nsim, iseed, init.value.theta_1, init.value.theta_2,init.value.theta_3) {
  per=0
  per2=0
  per3=0
  
  # Define empty tables to store the results of estimated CIF 
  # (nrow  = nr of iterations, ncol = nr of time points evaluated)
  # format: results.CiA.kj
  # i=1,2: C1 is competing risk 1 and C2 is competing risk 2
  # A= "","I","E","R", "C": with "" the two-step estimator, "I" the independent estimator,
  #    "E" the naive estimator, "C" a nonparametric estimator and "R" using the true parameter values
  # k=0,1: value of binary variable Z
  # j=0,1: value of binary variable W
  
  results.C1.11 = c() 
  results.C2.11 = c()
  results.C1I.11 = c()
  results.C2I.11 = c()
  results.C1E.11 = c()
  results.C2E.11 = c()
  results.C1R.11 = c()
  results.C2R.11 = c()
  results.C1C.11 = c()
  results.C2C.11 = c()
  
  results.C1.10 = c()
  results.C2.10 = c()
  results.C1I.10 = c()
  results.C2I.10 = c()
  results.C1E.10 = c()
  results.C2E.10 = c()
  results.C1R.10 = c()
  results.C2R.10 = c()
  results.C1C.10 = c()
  results.C2C.10 = c()
  
  results.C1.01 = c()
  results.C2.01 = c()
  results.C1I.01 = c()
  results.C2I.01 = c()
  results.C1E.01 = c()
  results.C2E.01 = c()
  results.C1R.01 = c()
  results.C2R.01 = c()
  results.C1C.01 = c()
  results.C2C.01 = c()
  
  results.C1.00 = c()
  results.C2.00 = c()
  results.C1I.00 = c()
  results.C2I.00 = c()
  results.C1E.00 = c()
  results.C2E.00 = c()
  results.C1R.00 = c()
  results.C2R.00 = c()
  results.C1C.00 = c()
  results.C2C.00 = c()
  
  
  for (i in 1:nsim) {
    # i = 1 # for testing
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,2,2) # generate data
    
    Y = data[,1]
    d1 = data[,2]
    d2 = data[,3]
    d3 = data[,4]
    da = data[,5]
    X = data[,(7:(parl+3))]
    Z = data[,parl+4]
    W = data[,parl+5]
    XandW = cbind(data[,6],X,W)
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    # Estimated V
    M = cbind(data[,6:(4+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,6:(5+parl)]
    
    # True value for V
    MrealV = cbind(data[,6:(4+parl)],data[,ncol(data)])
    
    per=per+table(d1)[2] # percentage T1 (C1)
    per2=per2+table(d2)[2] # percentage T2 (C2)
    per3=per3+table(d3)[2] # percentage dependent censoring
    
    # Model assuming independence between risks and independent censoring
    # (to generate starting values for the other models)
    
    # Assign starting values:
    # - beta1 = zero-vector
    # - beta2 = zero-vector
    # - beta3 = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - sigma3 = 1 
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    # - theta_3 = init.value.theta_3
    # - rho values are all assumed zero (no starting value necessary)
    
    init = c(rep(0,totparl), 1, 1,1, init.value.theta_1, init.value.theta_2, init.value.theta_3)
    parhat.init = nloptr(x0=c(init),eval_f=LikI.bis,Y=Y,d1=d1,d2=d2,d3=d3,da=da,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,1e-5, 0,0,0),ub=c(rep(Inf,totparl),Inf,Inf,Inf,2, 2,2),
                         eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
   
    # Independent model: Assuming independent censoring but no independence between risks

    # Assign starting values:
    # - beta1 = zero-vector
    # - beta2 = zero-vector
    # - beta3 = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - sigma3 = 1 
    # - rho12 = 0
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    # - theta_3 = init.value.theta_3
    
    # add starting value for rho12
    init<-c(parhat.init[1:(length(parhat.init)-3)],0,parhat.init[length(parhat.init)-2],parhat.init[length(parhat.init)-1],parhat.init[length(parhat.init)])
    
    parhat1 = nloptr(x0=c(init),eval_f=LikI.cmprsk,Y=Y,d1=d1,d2=d2,d3=d3,da=da,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,1e-5,-0.99, 0,0,0),ub=c(rep(Inf,totparl),Inf,Inf,Inf,0.99,2, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for V from data matrix
    ME = M[,-ncol(M)]
    
    # Remove coefficients for V in the vector parhat1. Add starting value for rho13 and rho23.
    # The final vector will be of the form:
    # [1:3] : beta1
    # [4:6] : beta2
    # [7:9] : beta3
    # [10]  : sigma1
    # [11]  : sigma2
    # [12]  : sigma3
    # [13]  : rho12
    # [14]  : rho13
    # [15]  : rho23
    # [16]  : theta_1
    # [17]  : theta_2
    # [18]  : theta_3
    
    # Remove coefficients for V 
    initE = parhat.init[-parl]
    initE = initE [-(2*parl-1)]
    initE = initE [-(3*parl-2)]
    
    # Append theta's to initE and replace the original theta_1 and theta_2 (now fourth and fifth-to-last
    # element) with the initial value for rho13 and rho23.
    initE = c(initE[-length(initE)],initE[length(initE)-2],initE[length(initE)-2],initE[length(initE)-1],initE[length(initE)])
    initE[length(initE) - 3] <- 0
    initE[length(initE) - 4] <- 0
    initE[length(initE) - 5] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF.cmprsk,Y=Y,d1=d1,d2=d2,d3=d3,da=da,M=ME,lb=c(rep(-Inf,(totparl-3)),1e-05,1e-5,1e-5,-0.99,-0.99,-0.99,0,0,0),ub=c(rep(Inf,(totparl-3)),Inf,Inf,Inf,0.99,0.99,0.99,2,2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=1000,"xtol_abs"=rep(1.0e-30)))$solution
    

    # Model with estimated V
    
    # Assign starting values
    # - beta1 (4 params) = First 4 params of parhat1
    # - beta2 (4 params) = Next 4 params of parhat1
    # - beta3 (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[13]
    # - sigma2 = parhat1[14]
    # - sigma3 = parhat1[15]
    # - rho12 = parhat1[16]
    # - rho13 = 0
    # - rho23 = 0
    # - theta_1 = parhat1[17]
    # - theta_2 = parhat1[18]
    # - theta_3 = parhat[19]
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-2],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 3] <- 0
    initd[length(initd) - 4] <- 0

    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhat = nloptr(x0=initd,eval_f=LikF.cmprsk,Y=Y,d1=d1,d2=d2,d3=d3,da=da,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,1e-5,-0.99,-0.99,-0.99,0,0,0),ub=c(rep(Inf,totparl),Inf,Inf,Inf,0.99,0.99,0.99,2,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    
    
    ####### Estimate CIF for different time-points (always assuming X=1) #######
    
    Time <- seq(from=1,to=200,by=0.5) # Times to evaluate
    
    # Z1-W1
    C1.11<-c() # create empty list to save the results for each time
    C2.11<-c()
    C1I.11<-c()
    C2I.11<-c()
    C1E.11<-c()
    C2E.11<-c()
    C1R.11<-c()
    C2R.11<-c()
    
    options(warn=-1)
    
    XandW.11<-c(1,1,1) #XandW matrix
    Z.11<-1 # Value for Z
    V.est.11 <- (1-Z.11)*((1+exp(XandW.11%*%gammaest))*log(1+exp(XandW.11%*%gammaest))-(XandW.11%*%gammaest)*exp(XandW.11%*%gammaest))-Z.11*((1+exp(-(XandW.11%*%gammaest)))*log(1+exp(-(XandW.11%*%gammaest)))+(XandW.11%*%gammaest)*exp(-(XandW.11%*%gammaest))) #Estimated control function
    V.true.11 <- (1-Z.11)*((1+exp(XandW.11%*%parN[[5]]))*log(1+exp(XandW.11%*%parN[[5]]))-(XandW.11%*%parN[[5]])*exp(XandW.11%*%parN[[5]]))-Z.11*((1+exp(-(XandW.11%*%parN[[5]])))*log(1+exp(-(XandW.11%*%parN[[5]])))+(XandW.11%*%parN[[5]])*exp(-(XandW.11%*%parN[[5]]))) #True control function
    M.11<-c(1,1,Z.11,V.est.11) # M matrix
    M.true.11 <- c(1,1,Z.11,V.true.11) # true M matrix
    for (i in 1:length(Time)){ #calculate CIF
      
      # Two-step estimator
      C1.11[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.11)
      C2.11[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.11)  
      
      # Independent estimator
      C1I.11[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.11)
      C2I.11[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.11)  
      
      # Naive estimator
      C1E.11[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,1))
      C2E.11[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,1))  
      
      # Using true parameter values
      C1R.11[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.11) 
      C2R.11[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.11) 
    }
    
    
    # Z1-W0
    C1.10<-c()
    C2.10<-c()
    C1I.10<-c()
    C2I.10<-c()
    C1E.10<-c()
    C2E.10<-c()
    C1R.10<-c()
    C2R.10<-c()
    
    XandW.10<-c(1,1,0)
    Z.10<-1
    V.est.10 <- (1-Z.10)*((1+exp(XandW.10%*%gammaest))*log(1+exp(XandW.10%*%gammaest))-(XandW.10%*%gammaest)*exp(XandW.10%*%gammaest))-Z.10*((1+exp(-(XandW.10%*%gammaest)))*log(1+exp(-(XandW.10%*%gammaest)))+(XandW.10%*%gammaest)*exp(-(XandW.10%*%gammaest)))
    V.true.10 <- (1-Z.10)*((1+exp(XandW.10%*%parN[[5]]))*log(1+exp(XandW.10%*%parN[[5]]))-(XandW.10%*%parN[[5]])*exp(XandW.10%*%parN[[5]]))-Z.10*((1+exp(-(XandW.10%*%parN[[5]])))*log(1+exp(-(XandW.10%*%parN[[5]])))+(XandW.10%*%parN[[5]])*exp(-(XandW.10%*%parN[[5]])))
    M.10<-c(1,1,Z.10,V.est.10)
    M.true.10 <- c(1,1,Z.10,V.true.10)
    for (i in 1:length(Time)){
      
      C1.10[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.10)
      C2.10[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.10)  
      
      C1I.10[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.10)
      C2I.10[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.10)  
      
      C1E.10[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,1))
      C2E.10[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,1))  
      
      C1R.10[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.10) 
      C2R.10[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.10) 
    }
    
    # Z0-W1
    C1.01<-c()
    C2.01<-c()
    C1I.01<-c()
    C2I.01<-c()
    C1E.01<-c()
    C2E.01<-c()
    C1R.01<-c()
    C2R.01<-c()
    
    XandW.01<-c(1,1,1)
    Z.01<-0
    V.est.01 <- (1-Z.01)*((1+exp(XandW.01%*%gammaest))*log(1+exp(XandW.01%*%gammaest))-(XandW.01%*%gammaest)*exp(XandW.01%*%gammaest))-Z.01*((1+exp(-(XandW.01%*%gammaest)))*log(1+exp(-(XandW.01%*%gammaest)))+(XandW.01%*%gammaest)*exp(-(XandW.01%*%gammaest)))
    V.true.01 <- (1-Z.01)*((1+exp(XandW.01%*%parN[[5]]))*log(1+exp(XandW.01%*%parN[[5]]))-(XandW.01%*%parN[[5]])*exp(XandW.01%*%parN[[5]]))-Z.01*((1+exp(-(XandW.01%*%parN[[5]])))*log(1+exp(-(XandW.01%*%parN[[5]])))+(XandW.01%*%parN[[5]])*exp(-(XandW.01%*%parN[[5]])))
    M.01<-c(1,1,Z.01,V.est.01)
    M.true.01 <- c(1,1,Z.01,V.true.01)
    for (i in 1:length(Time)){
      
      C1.01[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.01)
      C2.01[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.01)  
      
      C1I.01[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.01)
      C2I.01[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.01)  
      
      C1E.01[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,0))
      C2E.01[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,0))  
      
      C1R.01[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.01) 
      C2R.01[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.01) 
    }
    
    # Z0-W0
    C1.00<-c()
    C2.00<-c()
    C1I.00<-c()
    C2I.00<-c()
    C1E.00<-c()
    C2E.00<-c()
    C1R.00<-c()
    C2R.00<-c()
    
    XandW.00<-c(1,1,0)
    Z.00<-0
    V.est.00 <- (1-Z.00)*((1+exp(XandW.00%*%gammaest))*log(1+exp(XandW.00%*%gammaest))-(XandW.00%*%gammaest)*exp(XandW.00%*%gammaest))-Z.00*((1+exp(-(XandW.00%*%gammaest)))*log(1+exp(-(XandW.00%*%gammaest)))+(XandW.00%*%gammaest)*exp(-(XandW.00%*%gammaest)))
    V.true.00 <- (1-Z.00)*((1+exp(XandW.00%*%parN[[5]]))*log(1+exp(XandW.00%*%parN[[5]]))-(XandW.00%*%parN[[5]])*exp(XandW.00%*%parN[[5]]))-Z.00*((1+exp(-(XandW.00%*%parN[[5]])))*log(1+exp(-(XandW.00%*%parN[[5]])))+(XandW.00%*%parN[[5]])*exp(-(XandW.00%*%parN[[5]])))
    M.00<-c(1,1,Z.00,V.est.00)
    M.true.00 <- c(1,1,Z.00,V.true.00)
    for (i in 1:length(Time)){
      
      C1.00[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.00)
      C2.00[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[19],theta2=parhat[20],sigma1=parhat[13],sigma2=parhat[14],rho12=parhat[16],beta1=parhat[1:4],beta2=parhat[5:8],M=M.00)  
      
      C1I.00[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.00)
      C2I.00[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat1[17],theta2=parhat1[18],sigma1=parhat1[13],sigma2=parhat1[14],rho12=parhat1[16],beta1=parhat1[1:4],beta2=parhat1[5:8],M=M.00)  
      
      C1E.00[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,0))
      C2E.00[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[16],theta2=parhatE[17],sigma1=parhatE[10],sigma2=parhatE[11],rho12=parhatE[13],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,1,0))  
      
      C1R.00[i] <- integrate(integral1,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.00) 
      C2R.00[i] <- integrate(integral2,-Inf,log(Time[i]), theta1=parN[[4]][7],theta2=parN[[4]][8],sigma1=parN[[4]][1],sigma2=parN[[4]][2],rho12=parN[[4]][4],beta1=parN[[1]],beta2=parN[[2]],M=M.true.00) 
    }
    
    options(warn=0)
    
    # Bind the results from this iteration to tables
    
    results.C1.11 = rbind(results.C1.11,unlist(C1.11))
    results.C2.11 = rbind(results.C2.11,unlist(C2.11))
    results.C1I.11 = rbind(results.C1I.11,unlist(C1I.11))
    results.C2I.11 = rbind(results.C2I.11,unlist(C2I.11))
    results.C1E.11 = rbind(results.C1E.11,unlist(C1E.11))
    results.C2E.11 = rbind(results.C2E.11,unlist(C2E.11))
    results.C1R.11 = rbind(results.C1R.11,unlist(C1R.11))
    results.C2R.11 = rbind(results.C2R.11,unlist(C2R.11))
    
    results.C1.10 = rbind(results.C1.10,unlist(C1.10))
    results.C2.10 = rbind(results.C2.10,unlist(C2.10))
    results.C1I.10 = rbind(results.C1I.10,unlist(C1I.10))
    results.C2I.10 = rbind(results.C2I.10,unlist(C2I.10))
    results.C1E.10 = rbind(results.C1E.10,unlist(C1E.10))
    results.C2E.10 = rbind(results.C2E.10,unlist(C2E.10))
    results.C1R.10 = rbind(results.C1R.10,unlist(C1R.10))
    results.C2R.10 = rbind(results.C2R.10,unlist(C2R.10))
    
    results.C1.01 = rbind(results.C1.01,unlist(C1.01))
    results.C2.01 = rbind(results.C2.01,unlist(C2.01))
    results.C1I.01 = rbind(results.C1I.01,unlist(C1I.01))
    results.C2I.01 = rbind(results.C2I.01,unlist(C2I.01))
    results.C1E.01 = rbind(results.C1E.01,unlist(C1E.01))
    results.C2E.01 = rbind(results.C2E.01,unlist(C2E.01))
    results.C1R.01 = rbind(results.C1R.01,unlist(C1R.01))
    results.C2R.01 = rbind(results.C2R.01,unlist(C2R.01))
    
    results.C1.00 = rbind(results.C1.00,unlist(C1.00))
    results.C2.00 = rbind(results.C2.00,unlist(C2.00))
    results.C1I.00 = rbind(results.C1I.00,unlist(C1I.00))
    results.C2I.00 = rbind(results.C2I.00,unlist(C2I.00))
    results.C1E.00 = rbind(results.C1E.00,unlist(C1E.00))
    results.C2E.00 = rbind(results.C2E.00,unlist(C2E.00))
    results.C1R.00 = rbind(results.C1R.00,unlist(C1R.00))
    results.C2R.00 = rbind(results.C2R.00,unlist(C2R.00))
    
    
    
    
    # Nonparametric estimator of CIF
    data1 = data[data[,7]==1,] #value of X is equal to 1
    Y = exp(data1[,1])
    d1 = data1[,2]
    d2 = data1[,3]
    d3 = data1[,4]
    da = data1[,5]
    D = d1+2*d2 #1 of C1, 2 if C2 and 0 if censored (dependent or administrative)
    Z = data1[,parl+4]
    W = data1[,parl+5]
    group=factor(Z+2*W) #0 (Z=0,W=0), 1(Z=1,W=0), 2(Z=0,W=1), 3(Z=1,W=1)
      
    fit = timepoints(cuminc(Y,D,group,cencode=0),Time) # nonparametric estimation
      
    results.C1C.11 = rbind(results.C1C.11,fit$est[4,])
    results.C2C.11 = rbind(results.C2C.11,fit$est[8,])
    results.C1C.10 = rbind(results.C1C.10,fit$est[2,])
    results.C2C.10 = rbind(results.C2C.10,fit$est[6,])
    results.C1C.01 = rbind(results.C1C.01,fit$est[3,])
    results.C2C.01 = rbind(results.C2C.01,fit$est[7,])
    results.C1C.00 = rbind(results.C1C.00,fit$est[1,])
    results.C2C.00 = rbind(results.C2C.00,fit$est[5,])

  }
  # Our tables now contain the estimated CIF at different time points (columns)
  # for every iteration (rows)
  
  print(per/(n*nsim))     # percentage of C1
  print(per2/(n*nsim))    # percentage of C2
  print(per3/(n*nsim))    # percentage of dependent censoring
  
  # We will now create two tables:
  # results.mean contains the estimated mean CIF 
  #   - columns are time points
  #   - different row for every estimator, every competing risk and every Z-W combination
  # results.RMSE contains the estimated RMSE
  #   - columns are time points
  #   - different row for every estimator, every competing risk and every Z-W combination
  
  # Mean estimated CIF
  results.mean<-Time
  
  results.mean<-rbind(results.mean,apply(results.C1.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1I.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2I.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1E.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2E.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1R.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2R.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1C.11,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2C.11,2,mean))
  
  results.mean<-rbind(results.mean,apply(results.C1.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1I.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2I.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1E.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2E.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1R.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2R.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1C.10,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2C.10,2,mean))
  
  results.mean<-rbind(results.mean,apply(results.C1.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1I.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2I.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1E.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2E.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1R.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2R.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1C.01,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2C.01,2,mean))
  
  results.mean<-rbind(results.mean,apply(results.C1.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1I.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2I.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1E.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2E.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1R.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2R.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C1C.00,2,mean))
  results.mean<-rbind(results.mean,apply(results.C2C.00,2,mean))
  
  # Estimated RMSE
  results.RMSE<-Time
  
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1.11-results.C1R.11)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2.11-results.C2R.11)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1I.11-results.C1R.11)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2I.11-results.C2R.11)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1E.11-results.C1R.11)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2E.11-results.C2R.11)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1C.11-results.C1R.11)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2C.11-results.C2R.11)^2,2,mean)))
  
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1.10-results.C1R.10)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2.10-results.C2R.10)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1I.10-results.C1R.10)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2I.10-results.C2R.10)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1E.10-results.C1R.10)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2E.10-results.C2R.10)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1C.10-results.C1R.10)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2C.10-results.C2R.10)^2,2,mean)))
  
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1.01-results.C1R.01)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2.01-results.C2R.01)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1I.01-results.C1R.01)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2I.01-results.C2R.01)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1E.01-results.C1R.01)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2E.01-results.C2R.01)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1C.01-results.C1R.01)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2C.01-results.C2R.01)^2,2,mean)))
  
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1.00-results.C1R.00)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2.00-results.C2R.00)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1I.00-results.C1R.00)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2I.00-results.C2R.00)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1E.00-results.C1R.00)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2E.00-results.C2R.00)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C1C.00-results.C1R.00)^2,2,mean)))
  results.RMSE<-rbind(results.RMSE,sqrt(apply((results.C2C.00-results.C2R.00)^2,2,mean)))
  
  rownames(results.mean)<-c('time','C1.11','C2.11','C1I.11','C2I.11','C1E.11','C2E.11','C1R.11','C2R.11','C1C.11','C2C.11','C1.10','C2.10','C1I.10','C2I.10','C1E.10','C2E.10','C1R.10','C2R.10','C1C.10','C2C.10','C1.01','C2.01','C1I.01','C2I.01','C1E.01','C2E.01','C1R.01','C2R.01','C1C.01','C2C.01','C1.00','C2.00','C1I.00','C2I.00','C1E.00','C2E.00','C1R.00','C2R.00','C1C.00','C2C.00')
  rownames(results.RMSE)<-c('time','C1.11','C2.11','C1I.11','C2I.11','C1E.11','C2E.11','C1C.11','C2C.11','C1.10','C2.10','C1I.10','C2I.10','C1E.10','C2E.10','C1C.10','C2C.10','C1.01','C2.01','C1I.01','C2I.01','C1E.01','C2E.01','C1C.01','C2C.01','C1.00','C2.00','C1I.00','C2I.00','C1E.00','C2E.00','C1C.00','C2C.00')
  
  # Make nice Latex table
  xtab.mean = xtable(results.mean)
  
  # set to 3 significant digits
  digits(xtab.mean) = rep(3,400)
  
  header= c("mean")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab.mean,file=paste0("YJ_mean_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")

  # Make nice Latex table
  xtab.RMSE = xtable(results.RMSE)
  
  # set to 3 significant digits
  digits(xtab.RMSE) = rep(3,400)
  header= c("RMSE")

  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab.RMSE,file=paste0("YJ_RMSE_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  return(list(results.mean,results.RMSE)) # return tables containing mean and RMSE
  
}



########################## Data application cmprsk ##############################

# We can try to estimate the models gamma using
# 1. Firth's method 
# 2. putting V=0 for W=0 (control group) and estimating V for W=1 (screening group)
#    depending on the value of X (age) 


# Firth's method

DataApplication_cmprsk.Firth <- function(data, init.value.theta_1, init.value.theta_2) {
  
  
  n = nrow(data)
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  intercept = data[,4]
  X = data[,5:(parl + 1)]
  Z = data[,parl+2]
  W = data[,parl+3]
  XandW = cbind(intercept,X,W)

  # Estimate V using Firth regression
  gammaest <- summary(glm(as.factor(Z) ~ -1 + XandW, family = "binomial",method="brglmFit"))$coefficients
  gammaest <- gammaest[,1]

  ##############################################################################
  
  V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))

  # Create matrix of X, Z and V.
  M <- cbind(data[,4:(2+parl)],V)
  
  # Create matrix of X, Z and W.
  MnoV = data[,4:(3+parl)]
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  
  # Independent model for starting values sigmas and theta.
  # We use the likelihood from "Functions_ad" as we have only two variables (two competing risks)
  # The function LikI.cmprsk on the contrary assumes three variables (two competing risks + dependent censoring)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  print("independent")
  #
  # Our model: Taking into account Z is likely a confounded variable and that T1 
  #            and T2 are dependent
  #
  
  # Create vector of initial values for the estimation of the parameters using the
  # parhat1. The final vector will be of the following form. 
  # [1:4]  : beta1
  # [5:8]  : beta2
  # [9]    : sigma1
  # [10]   : sigma2
  # [11]   : rho
  # [12]   : theta_1
  # [13]   : theta_2
  
  # initial values: add starting value for rho
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  
  # We again use the function LikF and not LikF.cmprsk
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  print("two-step")
  
  parhatG = c(parhat,as.vector(gammaest))
  
  Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  
  # Select part of variance matrix pertaining to beta1, beta2, var1, var2, rho and theta
  # (i.e. H_delta).
  H = Hgamma[1:length(initd),1:length(initd)]
  HI = ginv(H)
  
  Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
  
  prodvec = XandW[,1]
  
  for (i in 1:parlgamma) {
    for (j in 2:parlgamma) {
      if (i<=j){
        prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
      }
    }
  }
  
  secder=t(-dlogis(XandW%*%gammaest))%*%prodvec
  
  WM = secder[1:parlgamma]
  for (i in 2:parlgamma) {
    newrow<-secder[c(i,(i+2):(i+parlgamma))]
    WM<-rbind(WM,newrow) 
  }
  
  WMI = ginv(WM)
  #Wh <- diag(as.vector(plogis(XandW)*(1-plogis(XandW))),nrow=n)
  #h <- Wh^(1/2)%*%XandW%*%ginv(t(XandW)%*%Wh%*%XandW)%*%t(XandW)%*%Wh^(1/2)
  #diffvec = Z-plogis(XandW%*%gammaest)+diag(h)*(1/2-plogis(XandW%*%gammaest))
  diffvec = Z-plogis(XandW%*%gammaest)
  
  mi = c()
  
  for(i in 1:n){
    newrow<-diffvec[i,]%*%XandW[i,]
    mi = rbind(mi,newrow)
  }
  
  psii = -WMI%*%t(mi)
  
  gi = c()
  
  for (i in 1:n)
  {
    J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    gi = rbind(gi,c(J1))
  }
  
  gi = t(gi)
  
  partvar = gi + Vargamma%*%psii
  
  Epartvar2 = (partvar%*%t(partvar))
  
  totvarex = HI%*%Epartvar2%*%t(HI)
  
  se = sqrt(abs(diag(totvarex)))
  
  # Delta method variance
  
  se_s1 = 1/parhat[totparl+1]*se[totparl+1]
  se_s2 = 1/parhat[totparl+2]*se[totparl+2]
  
  # Conf. interval for transf. sigma's
  
  st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
  st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
  
  # Back transform
  
  s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
  
  # Confidence interval for rho
  
  zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
  se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
  zt_l = zt-1.96*(se_z)
  zt_u = zt+1.96*(se_z)
  
  # Back transform
  
  r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
  r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
  
  # Confidence interval for theta
  
  rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
  rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
  rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
  rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
  
  # Matrix with all confidence intervals
  EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l, rtheta2_l),ncol=1),
              matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u,rtheta2_u), ncol=1))
  
  

  # Naive model: assuming Z is an unconfounded variable, but including dependence
  #              between two competing risks.
  #
  
  # remove data for v from data matrix
  ME = M[,-ncol(M)]
  
  # Remove coefficients for v in the vector parhat1. Add starting value for rho.
  # The final vector will be of the form:
  # [1:3]  : beta1
  # [5:6]  : beta2
  # [7]    : sigma1
  # [8]    : sigma2
  # [9]    : rho
  # [10]   : theta_1
  # [11]   : theta_2
  
  # Remove coefficients for v
  initE = parhat1[-parl]
  initE = initE[-(2*parl-1)]
  
  # Append theta to initE and replace the original theta (now second-to-last
  # element) with the initial value for rho.
  initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
  initE[length(initE) - 2] <- 0
  
  # Estimate the parameters
  # We again use the function LikF and not LikF.cmprsk
  parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,0.99,2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  print("naive")
 
  H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  H1I = ginv(H1)
  se1 = sqrt(abs(diag(H1I)));
  
  t_s1 = 1/parhatE[totparl-1]*se1[totparl-1]
  t_s2 = 1/parhatE[totparl]*se1[totparl]
  
  # Conf. interval for transf. sigma's
  
  ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
  ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
  
  # Back transform
  
  S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
  
  # Confidence interval for rho
  
  z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1])))     # Fisher's z transform
  se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
  z1t_l = z1t-1.96*(se1_z)
  z1t_u = z1t+1.96*(se1_z)
  
  # Back transform
  
  r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
  r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
  
  # Confidence interval for theta
  
  r1theta1_l <- parhatE[length(parhatE)-1] - 1.96 * se1[length(parhatE)-1]
  r1theta1_u <- parhatE[length(parhatE)-1] + 1.96 * se1[length(parhatE)-1]
  r1theta2_l <- parhatE[length(parhatE)] - 1.96 * se1[length(parhatE)]
  r1theta2_u <- parhatE[length(parhatE)] + 1.96 * se1[length(parhatE)]
  
  # Matrix of all the confidence intervals
  EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l, r1theta1_l, r1theta2_l),ncol=1),
              matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u, r1theta1_u, r1theta2_u),ncol=1)) 
  
  
  # Standard errors independent model
  
  parhatGI = c(parhat1,as.vector(gammaest))
  
  HgammaI = hessian(LikIGamma2,parhatGI,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  
  HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
  HIInd = ginv(HInd)
  
  VargammaI = HgammaI[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
  
  giI = c()
  
  for (i in 1:n) {
    J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    giI = rbind(giI,c(J1I))
  }
  
  giI = t(giI)
  
  partvarI = giI + VargammaI%*%psii
  
  Epartvar2I = (partvarI%*%t(partvarI))
  
  totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
  
  seI = sqrt(abs(diag(totvarexI)))
  
  # Delta method variance
  
  se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
  se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
  
  # Conf. interval for transf. sigma's
  
  st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
  st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
  
  # Back transform
  
  s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI)
  
  # Confidence interval for theta
  
  rItheta1_l <- parhat1[length(parhat1)-1] - 1.96 * seI[length(parhat1)-1]
  rItheta1_u <- parhat1[length(parhat1)-1] + 1.96 * seI[length(parhat1)-1]
  rItheta2_l <- parhat1[length(parhat1)] - 1.96 * seI[length(parhat1)]
  rItheta2_u <- parhat1[length(parhat1)] + 1.96 * seI[length(parhat1)]
  
  EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI,rItheta1_l, rItheta2_l),ncol=1),
              matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI,rItheta1_u, rItheta2_u), ncol=1))
  
  
  # Results of model assuming confounding and dependence between C1 and C2.
  pvalue <- 2*pmin((1-pnorm(parhat/se)),pnorm(parhat/se))
  significant <- ifelse(pvalue < 0.10,
                        ifelse(pvalue < 0.05,
                               ifelse(pvalue < 0.01, "**", "*"),"."), "")
  results.confound_dep <- cbind(parhat, se, pvalue, EC1)
  colnames(results.confound_dep) <- c("Estimate", "St.Dev.", "p", "CI.lb", "CI.ub")
  rownames(results.confound_dep) <- namescoef
  
  summary <- data.frame(round(results.confound_dep, 4))
  summary$sign <- significant
  summary <- summary[,c(1:3, 6, 4:5)]
  summary
  
  # Results of naive model
  pvalue.naive <- 2*pmin((1-pnorm(parhatE/se1)),pnorm(parhatE/se1))
  significant.naive <- ifelse(pvalue.naive < 0.10,
                              ifelse(pvalue.naive < 0.05,
                                     ifelse(pvalue.naive < 0.01, "**", "*"),"."), "")
  results.naive <- cbind(parhatE, se1, pvalue.naive, EC2)
  colnames(results.naive) <- c("Estimate", "St.Dev.", "pvalue", "CI.lb", "CI.ub")
  rownames(results.naive) <- namescoef[-c(parl, totparl)]
  
  summary1 <- data.frame(round(results.naive, 4))
  summary1$sign <- significant.naive
  summary1 <- summary1[,c(1:3, 6, 4:5)]
  summary1
  
  # Results of independence model
  pvalue.indep <- 2*pmin((1-pnorm(parhat1/seI)),pnorm(parhat1/seI))
  significant.indep <- ifelse(pvalue.indep < 0.10,
                              ifelse(pvalue.indep < 0.05,
                                     ifelse(pvalue.indep < 0.01, "**", "*"),"."), "")
  results.indep <- cbind(parhat1, seI, pvalue.indep, EC4)
  colnames(results.indep) <- c("Estimate", "St.Dev.", "pvalue", "CI.lb", "CI.ub")
  rownames(results.indep) <- namescoef[-(length(namescoef) - 2)]
  
  summary2 <- data.frame(round(results.indep, 4))
  summary2$sign <- significant.indep
  summary2 <- summary2[,c(1:3, 6, 4:5)]
  summary2
  
  ## Create LaTeX tables of results
  xtab = xtable(summary, digits = 3)
  header= c("sample size",n,"Results 2-step_Estimation with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_2-step_Estimation_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary1, digits = 3)
  header= c("sample size",n,"Results naive model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_naive_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary2, digits = 3)
  header= c("sample size",n,"Results independence model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  
  
  
  ########### Estimate CIF for different time-points (assuming unstandardized age of 50)
  Time <- seq(from=1,to=8000,by=1) # times to evaluate
  
  # Z1-W1 (people who actually participated in the screening)
  C1.11<-c() # using our model
  C2.11<-c()
  C1E.11<-c() # using model without endogeneity
  C2E.11<-c()
  print("start integral calculation")
  
  options(warn=-1)
  
  XandW.11<-c(1,-0.1360154,1) #XandW matrix
  Z.11<-1 # Value for Z
  V.est.11 <- (1-Z.11)*((1+exp(XandW.11%*%gammaest))*log(1+exp(XandW.11%*%gammaest))-(XandW.11%*%gammaest)*exp(XandW.11%*%gammaest))-Z.11*((1+exp(-(XandW.11%*%gammaest)))*log(1+exp(-(XandW.11%*%gammaest)))+(XandW.11%*%gammaest)*exp(-(XandW.11%*%gammaest))) #Estimated control function
  M.11<-c(1,-0.1360154,Z.11,V.est.11) # M matrix
  for (i in 1:length(Time)){ #calculate CIF
    
    C1.11[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.11)
    C2.11[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.11)  
    
    C1E.11[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,1))
    C2E.11[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,1))  
  }
  
  results.screened <- rbind(Time,C1.11,C2.11,C1E.11,C2E.11) # estimated CIF at different times for screened women
  

  # Z0-W1 (People who are selected for screening but placed themselves in control group)
  C1.01<-c()
  C2.01<-c()
  C1E.01<-c()
  C2E.01<-c()
  
  XandW.01<-c(1,-0.1360154,1)
  Z.01<-0
  V.est.01 <- (1-Z.01)*((1+exp(XandW.01%*%gammaest))*log(1+exp(XandW.01%*%gammaest))-(XandW.01%*%gammaest)*exp(XandW.01%*%gammaest))-Z.01*((1+exp(-(XandW.01%*%gammaest)))*log(1+exp(-(XandW.01%*%gammaest)))+(XandW.01%*%gammaest)*exp(-(XandW.01%*%gammaest)))
  M.01<-c(1,-0.1360154,Z.01,V.est.01)
  for (i in 1:length(Time)){
    
    C1.01[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.01)
    C2.01[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.01)  
   
    C1E.01[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))
    C2E.01[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))  
  }
  
  results.not_screened <- rbind(Time,C1.01,C2.01,C1E.01,C2E.01)
  
  # Z0-W0 (People in control group)
  C1.00<-c()
  C2.00<-c()
  C1E.00<-c()
  C2E.00<-c()

    XandW.00<-c(1,-0.1360154,0)
  Z.00<-0
  V.est.00 <- (1-Z.00)*((1+exp(XandW.00%*%gammaest))*log(1+exp(XandW.00%*%gammaest))-(XandW.00%*%gammaest)*exp(XandW.00%*%gammaest))-Z.00*((1+exp(-(XandW.00%*%gammaest)))*log(1+exp(-(XandW.00%*%gammaest)))+(XandW.00%*%gammaest)*exp(-(XandW.00%*%gammaest)))
  M.00<-c(1,-0.1360154,Z.00,V.est.00)
  for (i in 1:length(Time)){
    
    C1.00[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.00)
    C2.00[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.00)  
 
    C1E.00[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))
    C2E.00[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))  
   }
  
  results.control <- rbind(Time,C1.00,C2.00,C1E.00,C2E.00)
  
  options(warn=0)
 
  
  return(list(results.screened, results.not_screened, results.control)) # return all estimated CIF
  
  
  
}

# As people in the control group cannot participate in the screening, we can opt to put
# V=0 if W=0 and only estimate V if W=1

DataApplication_cmprsk <- function(data, init.value.theta_1, init.value.theta_2) {
  
  
  n = nrow(data)
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  intercept = data[,4]
  X = data[,5:(parl + 1)]
  Z = data[,parl+2]
  W = data[,parl+3]
  Xandintercept = cbind(intercept,X)
  
  # Estimate V if W=1
  
  gammaest <- summary(glm(as.factor(Z[W==1]) ~ -1 + Xandintercept[W==1,], family = "binomial"))$coefficients
  gammaest <- gammaest[,1]
  

  ##############################################################################
  
  # Construct V
  V <- rep(0,length(Z))
  V[W==1] <-(1-Z[W==1])*((1+exp(Xandintercept[W==1,]%*%gammaest))*log(1+exp(Xandintercept[W==1,]%*%gammaest))-(Xandintercept[W==1,]%*%gammaest)*exp(Xandintercept[W==1,]%*%gammaest))-Z[W==1]*((1+exp(-(Xandintercept[W==1,]%*%gammaest)))*log(1+exp(-(Xandintercept[W==1,]%*%gammaest)))+(Xandintercept[W==1,]%*%gammaest)*exp(-(Xandintercept[W==1,]%*%gammaest)))
  
  # Create matrix of X, Z and V.
  M <- cbind(data[,4:(2+parl)],V)
  
  # Create matrix of X, Z and W.
  MnoV = data[,4:(3+parl)]
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  
  # Independent model for starting values sigmas and theta.
  # We use the likelihood from "Functions_ad" as we have only two variables (two competing risks)
  # The function LikI.cmprsk on the contrary assumes three variables (two competing risks + dependent censoring)
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  print("independent")
  #
  # Our model: Taking into account Z is likely a confounded variable and that two 
  #            competing risks are dependent
  #
  
  # Create vector of initial values for the estimation of the parameters using the
  # parhat1. The final vector will be of the following form. 
  # [1:4]  : beta1
  # [5:8]  : beta2
  # [9]    : sigma1
  # [10]   : sigma2
  # [11]   : rho
  # [12]   : theta_1
  # [13]   : theta_2
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  # We again use the function LikF and not LikF.cmprsk
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  print("two-step")
  
  parhatG = c(parhat,as.vector(gammaest))
  
  Hgamma = hessian(LikFG.app,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  
  # Select part of variance matrix pertaining to beta1, beta2, var1, var2, rho and theta
  # (i.e. H_delta).
  H = Hgamma[1:length(initd),1:length(initd)]
  HI = ginv(H)
 
  Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
  
  prodvec = Xandintercept[W==1,1]
  
  for (i in 1:parlgamma) {
    for (j in 2:parlgamma) {
      if (i<=j){
        prodvec<-cbind(prodvec,diag(Xandintercept[W==1,i]%*%t(Xandintercept[W==1,j])))
      }
    }
  }
  
  secder=t(-dlogis(Xandintercept[W==1,]%*%gammaest))%*%prodvec
  
  WM = secder[1:parlgamma]
  for (i in 1:(parlgamma-1)) {
    newrow<-secder[c(i+1,(i+2):(i+parlgamma))]
    WM<-rbind(WM,newrow) 
  }
  
  WMI = ginv(WM)
  diffvec = Z[W==1]-plogis(Xandintercept[W==1,]%*%gammaest)
  
  mi = c()
  
  for(i in 1:sum(W==1)){
    newrow<-diffvec[i,]%*%(Xandintercept[W==1,])[i,]
    mi = rbind(mi,newrow)
  }
  
  psii = -WMI%*%t(mi)
  
  gi = c()
  
  for (i in 1:n)
  {
    J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    gi = rbind(gi,c(J1))
  }
  
  gi = t(gi)
  
  # Put psi equal to zero if V did not have to be estimated (W=0)
  psi <- matrix(0,nrow=dim(Vargamma)[1],ncol=n)
  sub <- Vargamma%*%psii
  
  j=1
  for (i in 1:n){
    if (W[i]==1){
      psi[,i] <- sub[,j]
      j=j+1
    }
  }
  
  partvar = gi + psi
  
  Epartvar2 = (partvar%*%t(partvar))
  
  totvarex = HI%*%Epartvar2%*%t(HI)
  
  se = sqrt(abs(diag(totvarex)))
  
  # Delta method variance
  
  se_s1 = 1/parhat[totparl+1]*se[totparl+1]
  se_s2 = 1/parhat[totparl+2]*se[totparl+2]
  
  # Conf. interval for transf. sigma's
  
  st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
  st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
  
  # Back transform
  
  s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
  
  # Confidence interval for rho
  
  zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
  se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
  zt_l = zt-1.96*(se_z)
  zt_u = zt+1.96*(se_z)
  
  # Back transform
  
  r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
  r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
  
  # Confidence interval for theta
  
  rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
  rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
  rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
  rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
  
  # Matrix with all confidence intervals
  EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l, rtheta2_l),ncol=1),
              matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u,rtheta2_u), ncol=1))
  
  
  
  # Naive model: assuming Z is an unconfounded variable, but including dependence
  #              between competing risks.
  #
  
  # remove data for v from data matrix
  ME = M[,-ncol(M)]
  
  # Remove coefficients for v in the vector parhat1. Add starting value for rho.
  # The final vector will be of the form:
  # [1:6]  : beta
  # [7:12] : eta
  # [13]   : sigma1
  # [14]   : sigma2
  # [15]   : rho
  # [16]   : theta_1
  # [17]   : theta_2
  
  # Remove coefficients for v
  initE = parhat1[-parl]
  initE = initE[-(2*parl-1)]
  
  # Append theta to initE and replace the original theta (now second-to-last
  # element) with the initial value for rho.
  initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
  initE[length(initE) - 2] <- 0
  
  # Estimate the parameters
  parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,0.99,2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  print("naive")
  
  H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  H1I = ginv(H1)
  se1 = sqrt(abs(diag(H1I)));
  
  t_s1 = 1/parhatE[totparl-1]*se1[totparl-1]
  t_s2 = 1/parhatE[totparl]*se1[totparl]
  
  # Conf. interval for transf. sigma's
  
  ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
  ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
  
  # Back transform
  
  S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
  
  # Confidence interval for rho
  
  z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1])))     # Fisher's z transform
  se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
  z1t_l = z1t-1.96*(se1_z)
  z1t_u = z1t+1.96*(se1_z)
  
  # Back transform
  
  r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
  r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
  
  # Confidence interval for theta
  
  r1theta1_l <- parhatE[length(parhatE)-1] - 1.96 * se1[length(parhatE)-1]
  r1theta1_u <- parhatE[length(parhatE)-1] + 1.96 * se1[length(parhatE)-1]
  r1theta2_l <- parhatE[length(parhatE)] - 1.96 * se1[length(parhatE)]
  r1theta2_u <- parhatE[length(parhatE)] + 1.96 * se1[length(parhatE)]
  
  # Matrix of all the confidence intervals
  EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l, r1theta1_l, r1theta2_l),ncol=1),
              matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u, r1theta1_u, r1theta2_u),ncol=1)) 
  
  
  # Standard errors for independent model
  
  parhatGI = c(parhat1,as.vector(gammaest))
  
  HgammaI = hessian(LikIGamma.app,parhatGI,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  
  HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
  HIInd = ginv(HInd)
  
  VargammaI = HgammaI[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
  
  giI = c()
  
  for (i in 1:n) {
    J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
    giI = rbind(giI,c(J1I))
  }
  
  giI = t(giI)
  
  psiI <- matrix(0,nrow=dim(VargammaI)[1],ncol=n)
  subI <- VargammaI%*%psii
  
  j=1
  for (i in 1:n){
    if (W[i]==1){
      psiI[,i] <- subI[,j]
      j=j+1
    }
  }
  
  partvarI = giI + psiI

  
  Epartvar2I = (partvarI%*%t(partvarI))
  
  totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
  
  seI = sqrt(abs(diag(totvarexI)))
  
  # Delta method variance
  
  se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
  se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
  
  # Conf. interval for transf. sigma's
  
  st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
  st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
  
  # Back transform
  
  s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI)
  
  # Confidence interval for theta
  
  rItheta1_l <- parhat1[length(parhat1)-1] - 1.96 * seI[length(parhat1)-1]
  rItheta1_u <- parhat1[length(parhat1)-1] + 1.96 * seI[length(parhat1)-1]
  rItheta2_l <- parhat1[length(parhat1)] - 1.96 * seI[length(parhat1)]
  rItheta2_u <- parhat1[length(parhat1)] + 1.96 * seI[length(parhat1)]
  
  EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI,rItheta1_l, rItheta2_l),ncol=1),
              matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI,rItheta1_u, rItheta2_u), ncol=1))
  
  
  # Results of model assuming confounding and dependence between competing risks.
  pvalue <- 2*pmin((1-pnorm(parhat/se)),pnorm(parhat/se))
  significant <- ifelse(pvalue < 0.10,
                        ifelse(pvalue < 0.05,
                               ifelse(pvalue < 0.01, "**", "*"),"."), "")
  results.confound_dep <- cbind(parhat, se, pvalue, EC1)
  colnames(results.confound_dep) <- c("Estimate", "St.Dev.", "p", "CI.lb", "CI.ub")
  rownames(results.confound_dep) <- namescoef
  
  summary <- data.frame(round(results.confound_dep, 4))
  summary$sign <- significant
  summary <- summary[,c(1:3, 6, 4:5)]
  summary
  
  # Results of naive model
  pvalue.naive <- 2*pmin((1-pnorm(parhatE/se1)),pnorm(parhatE/se1))
  significant.naive <- ifelse(pvalue.naive < 0.10,
                              ifelse(pvalue.naive < 0.05,
                                     ifelse(pvalue.naive < 0.01, "**", "*"),"."), "")
  results.naive <- cbind(parhatE, se1, pvalue.naive, EC2)
  colnames(results.naive) <- c("Estimate", "St.Dev.", "pvalue", "CI.lb", "CI.ub")
  rownames(results.naive) <- namescoef[-c(parl, totparl)]
  
  summary1 <- data.frame(round(results.naive, 4))
  summary1$sign <- significant.naive
  summary1 <- summary1[,c(1:3, 6, 4:5)]
  summary1
  
  # Results of independence model
  pvalue.indep <- 2*pmin((1-pnorm(parhat1/seI)),pnorm(parhat1/seI))
  significant.indep <- ifelse(pvalue.indep < 0.10,
                              ifelse(pvalue.indep < 0.05,
                                     ifelse(pvalue.indep < 0.01, "**", "*"),"."), "")
  results.indep <- cbind(parhat1, seI, pvalue.indep, EC4)
  colnames(results.indep) <- c("Estimate", "St.Dev.", "pvalue", "CI.lb", "CI.ub")
  rownames(results.indep) <- namescoef[-(length(namescoef) - 2)]
  
  summary2 <- data.frame(round(results.indep, 4))
  summary2$sign <- significant.indep
  summary2 <- summary2[,c(1:3, 6, 4:5)]
  summary2
  
  ## Create LaTeX tables of results
  xtab = xtable(summary, digits = 3)
  header= c("sample size",n,"Results 2-step_Estimation with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_2-step_Estimation_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary1, digits = 3)
  header= c("sample size",n,"Results naive model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_naive_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary2, digits = 3)
  header= c("sample size",n,"Results independence model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  
  
  
  ########### Estimate CIF for different time-points (assuming age of 50)
  Time <- seq(from=1,to=8000,by=1) # Times to evaluate
  
  # Z1-W1 (people who actually participated in the screening)
  C1.11<-c() # using our model
  C2.11<-c()
  C1E.11<-c() # using model without endogeneity
  C2E.11<-c()
  print("start integral calculation")
  
  options(warn=-1)
  
  X.11<-c(1,-0.1360154) #X matrix
  Z.11<-1 # Value for Z
  
  V.est.11 <- (1-Z.11)*((1+exp(X.11%*%gammaest))*log(1+exp(X.11%*%gammaest))-(X.11%*%gammaest)*exp(X.11%*%gammaest))-Z.11*((1+exp(-(X.11%*%gammaest)))*log(1+exp(-(X.11%*%gammaest)))+(X.11%*%gammaest)*exp(-(X.11%*%gammaest)))
  M.11<-c(1,-0.1360154,Z.11,V.est.11) # M matrix
  for (i in 1:length(Time)){ #calculate CIF
    
    C1.11[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.11)
    C2.11[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.11)  
    
    C1E.11[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,1))
    C2E.11[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,1))  
  }
  
  results.screened <- rbind(Time,C1.11,C2.11,C1E.11,C2E.11) # results for women who participated in screening
  
  
  # Z0-W1 (People who are selected for screening but place themselves in control group)
  C1.01<-c()
  C2.01<-c()
  C1E.01<-c()
  C2E.01<-c()
  
  X.01<-c(1,-0.1360154)
  Z.01<-0
  V.est.01 <- (1-Z.01)*((1+exp(X.01%*%gammaest))*log(1+exp(X.01%*%gammaest))-(X.01%*%gammaest)*exp(X.01%*%gammaest))-Z.01*((1+exp(-(X.01%*%gammaest)))*log(1+exp(-(X.01%*%gammaest)))+(X.01%*%gammaest)*exp(-(X.01%*%gammaest)))
  M.01<-c(1,-0.1360154,Z.01,V.est.01)
  for (i in 1:length(Time)){
    
    C1.01[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.01)
    C2.01[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.01)  
    
    C1E.01[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))
    C2E.01[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))  
  }
  
  results.not_screened <- rbind(Time,C1.01,C2.01,C1E.01,C2E.01)
  
  # Z0-W0 (People in control group)
  C1.00<-c()
  C2.00<-c()
  C1E.00<-c()
  C2E.00<-c()
  
  
  X.00<-c(1,-0.1360154)
  Z.00<-0
  V.est.00 <- 0
  M.00<-c(1,-0.1360154,Z.00,V.est.00)
  for (i in 1:length(Time)){
    
    C1.00[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.00)
    C2.00[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhat[12],theta2=parhat[13],sigma1=parhat[9],sigma2=parhat[10],rho12=parhat[11],beta1=parhat[1:4],beta2=parhat[5:8],M=M.00)  
    
    C1E.00[i] <-integrate(integral1,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))
    C2E.00[i] <-integrate(integral2,-Inf,log(Time[i]), theta1=parhatE[10],theta2=parhatE[11],sigma1=parhatE[7],sigma2=parhatE[8],rho12=parhatE[9],beta1=parhatE[1:3],beta2=parhatE[4:6],M=c(1,-0.1360154,0))  
  }
  
  results.control <- rbind(Time,C1.00,C2.00,C1E.00,C2E.00)
  
  options(warn=0)
  
  
  return(list(results.screened, results.not_screened, results.control)) # return all results
  
  
  
}


