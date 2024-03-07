######################### functions Yeo-Johnson ################################

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


###################### data simulating function ################################

dat.sim.reg = function(n,par,iseed,Zbin,Wbin){
  
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
  A = runif(n,0,8)
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

# likelihood gamma application
LikGamma_BC = function(par,Y,M){ 
  W=as.matrix(M)
  gamma= as.matrix(par)
  k=ncol(W)
  tilde_W = W[,k]
  Xandintercept = W[,-k]
  tot = (plogis(Xandintercept%*%gamma)^Y)*((1-plogis(Xandintercept%*%gamma))^(tilde_W-Y))
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

LikF = function(par,Y,Delta,Xi,M){ 
  M = as.matrix(M)
  k = ncol(M)
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k])
  eta = as.matrix(par[v:l])
  sigma1 = par[l+1]
  sigma2 = par[l+2]
  rho = par[l+3]
  theta_1 = par[l+4]
  theta_2 = par[l+5]
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(M%*%beta))/sigma1 # b_T
  z2 = ((transY.C-(rho*sigma2/sigma1)*transY.T)-(M%*%eta-rho*(sigma2/sigma1)*(M%*%beta)))/(sigma2*((1-rho^2)^0.5)) #  term within Phi for T
  z3 = (transY.C-(M%*%eta))/sigma2 # b_C
  z4 = ((transY.T-(rho*sigma1/sigma2)*transY.C)-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5) #  term within Phi for C
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4))*DtransY.C)^Xi*(pbinorm(q1=-z1,q2=-z3,cov12=rho))^(1-(Delta+Xi)) # likelihood
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}

# Needed for Hessian matrix when Z is continuous

# - Does assume endogeneity
#   - Z is continuous
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

LikFG1 = function(par,Y,Delta,Xi,M){ 
  M=as.matrix(M)
  k = ncol(M)-2 # remove Z and W
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
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=Z-XandW%*%gamma
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  # likelihood for estimated parameters theta(beta,eta,alpha, lambda, sigma,rho, theta) and gamma
  z1 = (transY.T-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((transY.C-rho*sigma2/sigma1*transY.T)-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (transY.C-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((transY.T-rho*sigma1/sigma2*transY.C)-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*(((1/sigma2)*dnorm(z3)*(1-pnorm(z4))*DtransY.C)^Xi)*((pbinorm(q1=-z1,q2=-z3,cov12=rho))^(1-(Delta+Xi))) # likelihood
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn) # returns likelihood (the same as likF but including gamma as parameter)
}

# needed for Hessian matrix when Z is binary

# - Does assume endogeneity
#   - Z is binary
# - Does assume dependent censoring
# - Uses Yeo-Johnson transformation

LikFG2 = function(par,Y,Delta,Xi,M){ 
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
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  
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

LikFG_BC= function(par,Y,Delta,Xi,M){ 
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
  XandW=as.matrix(cbind(X,W))
  
  Vest <-(W-Z)*((1+exp(X%*%gamma))*log(1+exp(X%*%gamma))-(X%*%gamma)*exp(X%*%gamma))-Z*((1+exp(-(X%*%gamma)))*log(1+exp(-(X%*%gamma)))+(X%*%gamma)*exp(-(X%*%gamma)))
 
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


################################
# Independent model assumption #
################################

# - Doesn't assume endogeneity
# - Doesn't assume dependent censoring
# - Uses Yeo-Johnson transformation

LikI = function(par,Y,Delta,Xi,M){ 
  M=as.matrix(M)
  k = ncol(M) 
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k]) # parameters related to T
  eta = as.matrix(par[v:l]) # parameters related to C
  sigma1 = par[l+1] # sigma of T
  sigma2 = par[l+2] # sigma of C
  theta_1 = par[l+3]
  theta_2 = par[l+4]
  
  transY.T=YJtrans(Y,theta_1)
  DtransY.T=DYJtrans(Y,theta_1)
  transY.C=YJtrans(Y,theta_2)
  DtransY.C=DYJtrans(Y,theta_2)
  
  z1 = (transY.T-(M%*%beta))/sigma1 # term within Phi for T (rho = 0)
  z2 = (transY.C-(M%*%eta))/sigma2 # term withing Phi for C (rho = 0)
  
  # likelihood (when rho = 0)
  # tot gives sample size number of rows
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2))*DtransY.T)^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1))*DtransY.C)^Xi)*((pbinorm(q1=-z1,q2=-z2,cov12=0))^(1-(Delta+Xi)))
  p1 = pmax(tot,1e-100) 
  Logn = sum(log(p1)); # calculate loglikelihood
  return(-Logn)
}

# needed for Hessian matrix

# - Does assume endogeneity
#   - Z is continuous
# - Doesn't assume dependent censoring
# - Uses Yeo-Johnson transformation

LikIGamma1 = function(par,Y,Delta,Xi,M){ 
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
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=Z-(XandW%*%gamma) #estimating the V 
  
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

# no dependent censoring (Z discrete)

# - Does assume endogeneity
#   - Z is binary
# - Doesn't assume dependent censoring
# - Uses Yeo-Johnson transformation

LikIGamma2 = function(par,Y,Delta,Xi,M){ 
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
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  
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



LikIGamma_BC = function(par,Y,Delta,Xi,M){ 
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
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest <-(W-Z)*((1+exp(X%*%gamma))*log(1+exp(X%*%gamma))-(X%*%gamma)*exp(X%*%gamma))-Z*((1+exp(-(X%*%gamma)))*log(1+exp(-(X%*%gamma)))+(X%*%gamma)*exp(-(X%*%gamma)))
  
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

# Some notes:
#   - How to set initial value for theta?
#   - In Deresa, the CI for theta is just constructed as
#         theta-hat +/- 1.96 se(theta-hat).
#     Shouldn't we impose CI(theta) \subset [0,2] (using, f.e. fisher-z transf)?
#
#     --> An appropriate transformation does not exist?
#
#   - I didn't change anything to the code to compute the variances on lines
#     ~ 671 - 730 (kan zijn dat de lijnnummers niet meer exact overeen komen
#     omdat, zoals nu, code/uitleg aan het bijschrijven ben). Not sure if 
#     something should change there when also including theta.
#   - If the code is ran with n = 100 matrices become singular. Needing so many
#     observations for an identified model seems like a downside.
#   - We kunnen de naam van deze functie altijd nog veranderen. Voor 't moment 
#     dacht ik dat het wel handig was om een andere naam te gebruiken dan in de 
#     code van Gilles of Deresa
#   - Wordt LikGamma1 ergens gebruikt? Zijn er functies die nergens gebruikt
#     worden?
################################################################################

SimulationCI11_SaraIlias = function(n, nsim, iseed, init.value.theta_1, init.value.theta_2) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    # i = 1 # for testing
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,1,1)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(2+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for v from data matrix
    ME = M[,-ncol(M)]
    
    # Remove coefficients for v in the vector parhat1. Add starting value for rho.
    # The final vector will be of the form:
    # [1:3] : beta
    # [4:6] : eta
    # [7]   : sigma1
    # [8]   : sigma2
    # [9]   : rho
    # [10]  : theta_1
    # [11]  : theta_2
    
    # Remove coefficients for v
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    
    # Append theta's to initE and replace the original theta_1 (now third-to-last
    # element) with the initial value for rho.
    initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
    initE[length(initE) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,0.99,2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    # --> Take the log of the estimates, construct CI for the log-estimates and
    #     backtransform to the original estimates by exponentiating.
    
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
    r1_u = min(((exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    r1theta1_l <- parhatE[length(parhatE)-1] - 1.96 * se1[length(parhatE)-1]
    r1theta1_u <- parhatE[length(parhatE)-1] + 1.96 * se1[length(parhatE)-1]
    r1theta2_l <- parhatE[length(parhatE)] - 1.96 * se1[length(parhatE)]
    r1theta2_u <- parhatE[length(parhatE)] + 1.96 * se1[length(parhatE)]
    
    # Matrix of all the confidence intervals
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l, r1theta1_l, r1theta2_l),ncol=1),
                matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u, r1theta1_u, r1theta2_u),ncol=1)) 
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
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
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta's
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
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    # M-matrix: second derivative of m(W,Z,gamma)
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    # Inverse of M-matrix
    WMI = ginv(WM)
    
    # First derivative of m(W,Z,gamma)
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    # psi_i-matrix
    psii = -WMI%*%mi
    
    # h_l(S_i, gamma, delta)
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    # h_l(S, gamma, delta) + H_gamma %*% Psi_i
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
    r_u = min(((exp(2*zt_u)-1)/(exp(2*zt_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rtheta1_l <- parhat[length(parhat)-1] - 1.96 * se[length(parhat)-1]
    rtheta1_u <- parhat[length(parhat)-1] + 1.96 * se[length(parhat)-1]
    rtheta2_l <- parhat[length(parhat)] - 1.96 * se[length(parhat)]
    rtheta2_u <- parhat[length(parhat)] + 1.96 * se[length(parhat)]
    
    # Matrix with all confidence intervals
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = min(((exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)),1,na.rm=TRUE)
    
    # Confidence interval for theta
    
    rretheta1_l <- parhatre[length(parhatre)-1] - 1.96 * sere[length(parhatre)-1]
    rretheta1_u <- parhatre[length(parhatre)-1] + 1.96 * sere[length(parhatre)-1]
    rretheta2_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
    rretheta2_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta1_l, rretheta2_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta1_u, rretheta2_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    # We construct the vector with
    # [1:4]   : params for beta
    # [5:8]   : params for eta
    # [9]     : param for sigma1
    # [10]    : param for sigma2
    # [11]    : param for theta_1
    # [12]    : param for theta_2
    # [13,15] : params for (gamma_0, gamma_X, gamma_W)
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma1,parhatGI,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = HgammaI[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
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
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim)) 
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with no V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  
  # Remove parameters pertaining to V
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  # par0:
  # - [1:3] : beta
  # - [4:6] : eta
  # - [7]   : sigma1
  # - [8]   : sigma2
  # - [9]   : rho
  # - [10]  : theta_1
  # - [11]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results1[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+3)],2,sd)
  RMSE = sqrt(apply((results1[,1:(totparl+3)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results1[,((totparl+3) + 1):(2*(totparl+3))],2,mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,(totparl+3))
  datacp = results1[,(2*(totparl+3)+1):(4*(totparl+3))]
  for(i in 1:(totparl+3)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with real V
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results2[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results2[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results2[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Results of model with estimated V but independence
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2], parN[[3]][4], parN[[3]][5])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : theta_1
  # - [12]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates.
  Bias = apply(results3[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results3[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results3[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  CP = rep(0,totparl+4)
  datacp = results3[,(2*(totparl+4) + 1):(4*(totparl+4))]
  for (i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  
  # Make nice Latex table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,6)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab,file=paste0("YJ_estV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("YJ_noV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("YJ_realV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = c(namescoef[1:(length(namescoef)-3)],namescoef[length(namescoef)-1],namescoef[length(namescoef)])
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("YJ_IndEstV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}

# Some notes:
#   - Same notes as for SimulationCI11_SaraIlias
#   - The only difference between this code and Simulation CI11_SaraIlias is 
#     on line 1170 where data is generated (of course the names of the txt-files
#     storing the results is also different.)
################################################################################

SimulationCI12_SaraIlias = function(n, nsim, iseed, init.value.theta_1, init.value.theta_2) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    # i = 1 # for testing
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,1,2)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(2+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for v from data matrix
    ME = M[,-ncol(M)]
    
    # Remove coefficients for v in the vector parhat1. Add starting value for rho.
    # The final vector will be of the form:
    # [1:3] : beta
    # [4:6] : eta
    # [7]   : sigma1
    # [8]   : sigma2
    # [9]   : rho
    # [10]  : theta_1
    # [11]  : theta_2
    
    # Remove coefficients for v
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    
    # Append theta's to initE and replace the original theta_1 (now third-to-last
    # element) with the initial value for rho.
    initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
    initE[length(initE) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,0.99,2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    # --> Take the log of the estimates, construct CI for the log-estimates and
    #     backtransform to the original estimates by exponentiating.
    
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
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
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
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
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
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    # M-matrix: second derivative of m(W,Z,gamma)
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    # Inverse of M-matrix
    WMI = ginv(WM)
    
    # First derivative of m(W,Z,gamma)
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    # psi_i-matrix
    psii = -WMI%*%mi
    
    # h_l(S_i, gamma, delta)
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    # h_l(S, gamma, delta) + H_gamma %*% Psi_i
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
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    # Confidence interval for theta
    
    rretheta1_l <- parhatre[length(parhatre)-1] - 1.96 * sere[length(parhatre)-1]
    rretheta1_u <- parhatre[length(parhatre)-1] + 1.96 * sere[length(parhatre)-1]
    rretheta2_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
    rretheta2_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta1_l, rretheta2_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta1_u, rretheta2_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    # We construct the vector with
    # [1:4]   : params for beta
    # [5:8]   : params for eta
    # [9]     : param for sigma1
    # [10]    : param for sigma2
    # [11]    : param for theta_1
    # [12]    : param for theta_2
    # [13,15] : params for (gamma_0, gamma_X, gamma_W)
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma1,parhatGI,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = HgammaI[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
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
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim)) 
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with no V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  
  # Remove parameters pertaining to V
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  # par0:
  # - [1:3] : beta
  # - [4:6] : eta
  # - [7]   : sigma1
  # - [8]   : sigma2
  # - [9]   : rho
  # - [10]  : theta_1
  # - [11]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results1[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+3)],2,sd)
  RMSE = sqrt(apply((results1[,1:(totparl+3)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results1[,((totparl+3) + 1):(2*(totparl+3))],2,mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,(totparl+3))
  datacp = results1[,(2*(totparl+3)+1):(4*(totparl+3))]
  for(i in 1:(totparl+3)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with real V
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [10]  : theta_1
  # - [11]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results2[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results2[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results2[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Results of model with estimated V but independence
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2], parN[[3]][4], parN[[3]][5])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : theta
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates.
  Bias = apply(results3[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results3[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results3[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  CP = rep(0,totparl+4)
  datacp = results3[,(2*(totparl+4) + 1):(4*(totparl+4))]
  for (i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  
  # Make nice Latex table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,6)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab,file=paste0("YJ_estV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("YJ_noV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("YJ_realV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = c(namescoef[1:(length(namescoef)-3)],namescoef[length(namescoef)-1],namescoef[length(namescoef)])
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("YJ_IndEstV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}



# Some notes:
#   - Same notes as for SimulationCI11
#   - Only differences with SimulationCI11 are:
#       x The parameters in the data-generation line of code in beginning
#       x The way gammaest and V are constructed not much later
#       x When estimating Hgamma under parhatG, we need to use LikFG2 instead of 
#         LikFG1.
#       x The way in which (sum)secder is made a little further on also changes.
#       x And hence also WM, diffvec, psii, ... on the next lines
#       x When estimating HgammaI under parhatGI, we need to use LikGamma2
#         instead of LikGamma1.
#       x Change names of txt-files
################################################################################


SimulationCI21_SaraIlias = function(n, nsim, iseed, init.value.theta_1, init.value.theta_2) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    # i = 1 # for testing
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,2,1)
    
    Y = data[,1]
    Delta = data[,2]
    Xi=data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(2+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta = init.value.theta
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for v from data matrix
    ME = M[,-ncol(M)]
    
    # Remove coefficients for v in the vector parhat1. Add starting value for rho.
    # The final vector will be of the form:
    # [1:3] : beta
    # [4:6] : eta
    # [7]   : sigma1
    # [8]   : sigma2
    # [9]   : rho
    # [10]  : theta_1
    # [11]  : theta_2
    
    # Remove coefficients for v
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    
    # Append theta's to initE and replace the original theta_1 (now third-to-last
    # element) with the initial value for rho.
    initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
    initE[length(initE) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,0.99,2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    # --> Take the log of the estimates, construct CI for the log-estimates and
    #     backtransform to the original estimates by exponentiating.
    
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
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
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
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
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
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    # Confidence interval for theta
    
    rretheta1_l <- parhatre[length(parhatre)-1] - 1.96 * sere[length(parhatre)-1]
    rretheta1_u <- parhatre[length(parhatre)-1] + 1.96 * sere[length(parhatre)-1]
    rretheta2_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
    rretheta2_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta1_l, rretheta2_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta1_u, rretheta2_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    # We construct the vector with
    # [1:4]   : params for beta
    # [5:8]   : params for eta
    # [9]     : param for sigma1
    # [10]    : param for sigma2
    # [11]    : param for theta_1
    # [12]    : param for theta_2
    # [13,15] : params for (gamma_0, gamma_X, gamma_W)
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
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with no V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  
  # Remove parameters pertaining to V
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  # par0:
  # - [1:3] : beta
  # - [4:6] : eta
  # - [7]   : sigma1
  # - [8]   : sigma2
  # - [9]   : rho
  # - [10]  : theta_1
  # - [11]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results1[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+3)],2,sd)
  RMSE = sqrt(apply((results1[,1:(totparl+3)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results1[,((totparl+3) + 1):(2*(totparl+3))],2,mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,(totparl+3))
  datacp = results1[,(2*(totparl+3)+1):(4*(totparl+3))]
  for(i in 1:(totparl+3)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with real V
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results2[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results2[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results2[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Results of model with estimated V but independence
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2], parN[[3]][4], parN[[3]][5])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : theta_1
  # - [12]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates.
  Bias = apply(results3[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results3[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results3[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  CP = rep(0,totparl+4)
  datacp = results3[,(2*(totparl+4) + 1):(4*(totparl+4))]
  for (i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  
  # Make nice Latex table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,6)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab,file=paste0("YJ_estV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("YJ_noV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("YJ_realV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = c(namescoef[1:(length(namescoef)-3)],namescoef[length(namescoef)-1],namescoef[length(namescoef)])
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("YJ_IndEstV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}


# Some notes:
#   - Same notes as for SimulationCI21
#   - Only differences with SimulationCI21 are:
#       - data generation arguments
#       - Names of files
################################################################################


SimulationCI22_SaraIlias = function(n, nsim, iseed, init.value.theta_1, init.value.theta_2) {
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  per2=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim) {
    # i = 1 # for testing
    
    if (round(i %% (nsim/10)) == 0) {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,2,2)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(2+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V --> remove data for v from data matrix
    ME = M[,-ncol(M)]
    
    # Remove coefficients for v in the vector parhat1. Add starting value for rho.
    # The final vector will be of the form:
    # [1:3] : beta
    # [4:6] : eta
    # [7]   : sigma1
    # [8]   : sigma2
    # [9]   : rho
    # [10]  : theta_1
    # [11]  : theta_2
    
    # Remove coefficients for v
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    
    # Append theta's to initE and replace the original theta_1 (now third-to-last
    # element) with the initial value for rho.
    initE = c(initE[-length(initE)],initE[length(initE)-1],initE[length(initE)])
    initE[length(initE) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,(totparl-2)),Inf,Inf,0.99,2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,Xi=Xi,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    # --> Take the log of the estimates, construct CI for the log-estimates and
    #     backtransform to the original estimates by exponentiating.
    
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
    
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
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
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
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
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    # Confidence interval for theta
    
    rretheta1_l <- parhatre[length(parhatre)-1] - 1.96 * sere[length(parhatre)-1]
    rretheta1_u <- parhatre[length(parhatre)-1] + 1.96 * sere[length(parhatre)-1]
    rretheta2_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
    rretheta2_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta1_l, rretheta2_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta1_u, rretheta2_u), ncol=1))
    # Model with estimated V but assuming independence
    
    # We construct the vector with
    # [1:4]   : params for beta
    # [5:8]   : params for eta
    # [9]     : param for sigma1
    # [10]    : param for sigma2
    # [11]    : param for theta_1
    # [12]    : param for theta_2
    # [13:15] : params for (gamma_0, gamma_X, gamma_W)
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
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
  }
  
  print(per/(n*nsim))     #percentage of censoring
  print(per2/(n*nsim))
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with no V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  
  # Remove parameters pertaining to V
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  # par0:
  # - [1:3] : beta
  # - [4:6] : eta
  # - [7]   : sigma1
  # - [8]   : sigma2
  # - [9]   : rho
  # - [10]  : theta_1
  # - [11]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results1[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+3)],2,sd)
  RMSE = sqrt(apply((results1[,1:(totparl+3)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results1[,((totparl+3) + 1):(2*(totparl+3))],2,mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,(totparl+3))
  datacp = results1[,(2*(totparl+3)+1):(4*(totparl+3))]
  for(i in 1:(totparl+3)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with real V
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results2[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results2[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results2[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Results of model with estimated V but independence
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2], parN[[3]][4], parN[[3]][5])
  par0m = matrix(par0,nsim,(totparl+4),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : theta_1
  # - [12]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates.
  Bias = apply(results3[,1:(totparl+4)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+4)],2,sd)
  RMSE = sqrt(apply((results3[,1:(totparl+4)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results3[,((totparl+4)+1):(2*(totparl+4))],2, mean)
  
  CP = rep(0,totparl+4)
  datacp = results3[,(2*(totparl+4) + 1):(4*(totparl+4))]
  for (i in 1:(totparl+4)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  
  # Make nice Latex table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,6)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab,file=paste0("YJ_estV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("YJ_noV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("YJ_realV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = c(namescoef[1:(length(namescoef)-3)],namescoef[length(namescoef)-1],namescoef[length(namescoef)])
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("YJ_IndEstV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}


#################### SimulationCI11_SaraIlias_Simplified #######################

# A simplified version of the above function, only fitting the two-step 
# estimator model and oracle model. There is also the added functionality of 
# being able to split up the simulation in a specified amount of parts so that
# we do not have to run this lengthy simulation in one session.
SimulationCI11_SaraIlias_Simplified = function(n, nsim, iseed, 
                                               init.value.theta_1,
                                               init.value.theta_2,
                                               part.to.evaluate,
                                               number.of.parts) {
  
  # Some input validation
  if (nsim %% number.of.parts != 0) {
    stop("nsim needs to be a multiple of number.of.parts.")
  }
  if ((part.to.evaluate > number.of.parts) || (part.to.evaluate <= 0)) {
    stop("part.to.evaluate is not valid.")
  }
  
  # Create the appropriate set of i's to check
  part.size <- nsim / number.of.parts
  i.to.check <- 1:nsim
  i.to.check <- i.to.check[(part.size*(part.to.evaluate-1) + 1):(part.size*part.to.evaluate)]
  
  sum.estV = c()
  sum.realV = c()
  per = 0
  per2 = 0
  results.estV = c()
  results.realV = c()
  
  for (i in i.to.check) {
    
    cat((i - i.to.check[1] + 1)*100/length(i.to.check), "% Completion \n")
    
    data = dat.sim.reg(n,parN,iseed+i,1,1)
    
    Y = data[,1]
    Delta = data[,2]
    Xi = data[,3]
    X = data[,(5:(parl+1))]
    Z = data[,parl+2]
    W = data[,parl+3]
    XandW = cbind(data[,4],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    # Estimated V
    M = cbind(data[,4:(2+parl)],V)
    
    # No V (using W instead)
    MnoV = data[,4:(3+parl)]
    
    # True value for V
    MrealV = cbind(data[,4:(2+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[2]
    per2=per2+table(Xi)[2]
    
    # Assign starting values:
    # - beta = zero-vector
    # - eta = zero-vector
    # - sigma1 = 1
    # - sigma2 = 1 
    # - theta_1 = init.value.theta_1
    # - theta_2 = init.value.theta_2
    init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
    
    # Independent model for starting values sigma and theta.
    #
    # Note the difference with the version of Gilles: the likelihood function now
    # takes an extra argument (= theta), so the vector of initial values needs
    # to take this into account. Also the vectors for the lower -and upper bound
    # of the parameters ('lb' and 'ub') should take this into account. Note that
    # theta is a value between 0 and 2.
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with estimated V
    
    # Assign starting values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
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
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta's
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
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    # M-matrix: second derivative of m(W,Z,gamma)
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    # Inverse of M-matrix
    WMI = ginv(WM)
    
    # First derivative of m(W,Z,gamma)
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    # psi_i-matrix
    psii = -WMI%*%mi
    
    # h_l(S_i, gamma, delta)
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],Xi=Xi[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    # h_l(S, gamma, delta) + H_gamma %*% Psi_i
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
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l,rtheta1_l,rtheta2_l),ncol=1),
                matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u,rtheta1_u, rtheta2_u), ncol=1))
    
    # Model with real V
    
    # Retake vector with initial values
    # - beta (4 params) = First 4 params of parhat1
    # - eta (4 params) = Next 4 params of parhat1
    # - sigma1 = parhat1[9]
    # - sigma2 = parhat1[10]
    # - rho = 0
    # - theta_1 = parhat1[11]
    # - theta_2 = parhat1[12]
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    # Again we make sure to properly adapt the upper -and lower bound values of
    # theta.
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,Xi=Xi,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    # Confidence interval for theta
    
    rretheta1_l <- parhatre[length(parhatre)-1] - 1.96 * sere[length(parhatre)-1]
    rretheta1_u <- parhatre[length(parhatre)-1] + 1.96 * sere[length(parhatre)-1]
    rretheta2_l <- parhatre[length(parhatre)] - 1.96 * sere[length(parhatre)]
    rretheta2_u <- parhatre[length(parhatre)] + 1.96 * sere[length(parhatre)]
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l,rretheta1_l, rretheta2_l),ncol=1),
                matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u,rretheta1_u, rretheta2_u), ncol=1))
    
    results.estV = rbind(results.estV,c(parhat,se,c(t(EC1))))
    results.realV = rbind(results.realV,c(parhatre,sere,c(t(EC3))))
  }
  
  # For each simulation, we also record the percentage of censored and admin-
  # istratively censored observations. Note that the global percentages can be 
  # easily computed by taking the average over all of the separate percentages,
  # as, by design, n*length(i.to.check) will be the same no matter the value of
  # part.to.evaluate.
  percentage1 <- per/(n*length(i.to.check))     #percentage of censoring
  percentage2 <- per2/(n*length(i.to.check))
  
  df.estV <- data.frame(results.estV, row.names = i.to.check)
  df.realV <- data.frame(results.realV, row.names = i.to.check)
  df.percentage <- data.frame(per1 = percentage1, per2 = percentage2)
  
  write.csv(df.estV, file = paste0("Sim data/df_estV_", part.to.evaluate,
                                   "_out_of_", number.of.parts, ".csv"),
            row.names = FALSE)
  
  write.csv(df.realV, file = paste0("Sim data/df_realV_", part.to.evaluate,
                                    "_out_of_", number.of.parts, ".csv"),
            row.names = FALSE)
  
  write.csv(df.percentage, file = paste0("Sim data/df_percentage_", part.to.evaluate,
                                         "_out_of_", number.of.parts, ".csv"),
            row.names = FALSE)
  
}

########################### Summarize_results ##################################

# This function collects all data files from the above simulation and performs
# the final analysis on them.
Summarize_results = function() {
  
  #
  # Create full data sets
  #
  
  # Get all file names in 'Chess data' folder
  files <- list.files("Sim data")
  
  # Read all files starting with "df_estV_". Store them in a list object. Do the
  # same for files starting with "df_realV_" and "df_percentage_"
  data_estV <- list()
  data_realV <- list()
  data_percentage <- list()
  
  for (file_name in files) {
    if (grepl("df_estV_", file_name)) {
      data_estV[[length(data_estV) + 1]] <- read.csv(paste0("Sim data/", file_name))
    }
    if (grepl("df_realV_", file_name)) {
      data_realV[[length(data_realV) + 1]] <- read.csv(paste0("Sim data/", file_name))
    }
    if (grepl("df_percentage_", file_name)) {
      data_percentage[[length(data_percentage) + 1]] <- read.csv(paste0("Sim data/", file_name))
    }
  }
  
  # Create empty data frames
  results.estV <- data_estV[[1]]
  results.realV <- data_realV[[1]]
  results.percentage <- data_percentage[[1]]
  
  # Append all separate data frames
  for (i in 2:length(data_estV)) {
    results.estV <- rbind(results.estV, data_estV[[i]])
    results.realV <- rbind(results.realV, data_realV[[i]])
    results.percentage <- rbind(results.percentage, data_percentage[[i]])
  }
  
  #
  # Results of model with estimated V
  #
  
  # Put all parameters (except gamma) into a vector
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results.estV[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results.estV[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results.estV[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the parameter standard deviations
  MSD  = apply(results.estV[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results.estV[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  #
  # Model with real V
  #
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+5),byrow=TRUE)
  # par0:
  # - [1:4] : beta
  # - [5:8] : eta
  # - [9]   : sigma1
  # - [10]  : sigma2
  # - [11]  : rho
  # - [12]  : theta_1
  # - [13]  : theta_2
  #
  # - totparl = 8
  
  # Statistics on the parameter estimates
  Bias = apply(results.realV[,1:(totparl+5)]-par0m,2,mean)
  ESE = apply(results.realV[,1:(totparl+5)],2,sd)
  RMSE = sqrt(apply((results.realV[,1:(totparl+5)]-par0m)^2,2,mean))
  
  # Statistics on the standard deviation estimates
  MSD  = apply(results.realV[,((totparl+5)+1):(2*(totparl+5))],2, mean)
  
  # Statistics on the parameter CI's: for each parameter, check how many times the
  # true value is contained in the estimated confidence interval. We divide by
  # nsim to obtain a percentage.
  CP = rep(0,totparl+5)
  datacp = results.realV[,(2*(totparl+5)+1):(4*(totparl+5))]
  for(i in 1:(totparl+5)) {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum2 = summary2
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  
  # Make nice Latex table
  xtab = xtable(sum)
  
  # set to 3 significant digits
  digits(xtab) = rep(3,6)
  
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  # Save table code in .txt-file. Also add header row.
  print.xtable(xtab,file=paste0("YJ_estV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("YJ_realV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
}

########################### Data application JTPA ##############################

DataApplicationJPTA <- function(data, init.value.theta_1, init.value.theta_2,
                                multiple.starting.points = FALSE) {
  
  n = nrow(data)
  Y = data[, 1]
  Delta = data[, 2]
  Xi = data[, 3]
  X = data[,(5:(parl + 1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  XandW = as.matrix(cbind(data[,4],X,W))
  
  # Estimate V
  gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  
  # Create matrix of X, Z and V.
  M <- cbind(data[,4:(2+parl)],V)
  
  # Create matrix of X, Z and W.
  MnoV = data[,4:(3+parl)]
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  
  # Independent model for starting values sigmas and theta.
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  #
  # Our model: Taking into account Z is likely a confounded variable and that T 
  #            and C are dependent
  #
  
  # Create vector of initial values for the estimation of the parameters using the
  # parhat1. The final vector will be of the following form. 
  # [1:7]  : beta
  # [8:14] : eta
  # [15]   : sigma1
  # [16]   : sigma2
  # [17]   : rho
  # [18]   : theta_1
  # [19]   : theta_2
  
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
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
  }
  
  # Concatenate estimated coefficients vector for the models for transformed event
  # and censoring time with the vector of coefficients for the model
  # Z = \gamma * W + V. This results in the following vector.
  # [1:7]   : beta
  # [8:14]  : eta
  # [15]    : sigma1
  # [16]    : sigma2
  # [17]    : rho
  # [18]    : theta_1
  # [19]   : theta_2
  # [20:25] : gamma
  parhatG = c(parhat,as.vector(gammaest))
  
  Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=10, v=2, show.details=FALSE)) 
  
  # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
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
  WM<-rbind(WM,secder[c(2,5,6,7)])
  WM<-rbind(WM,secder[c(3,6,8,9)])
  WM<-rbind(WM,secder[c(4,7,9,10)])

  WMI = ginv(WM)
  
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
  
  
  #
  # Naive model: assuming Z is an unconfounded variable, but including dependence
  #              between T and C.
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
  
  #
  # Independence model: Model taking into account that Z is likely a confounded
  #                     variable but assuming independence between T and C.
  #
  
  # We construct the vector with
  # [1:7]   : beta
  # [8:14]  : eta
  # [15]    : sigma1
  # [16]    : sigma2
  # [17]    : theta_1
  # [18]    : theta_2
  # [19:24] : gamma
  
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
  
  
  # Results of model assuming confounding and dependence between T and C.
  pvalue <- 2*pmin((1-pnorm(parhat/se)),pnorm(parhat/se))
  significant <- ifelse(pvalue < 0.10,
                        ifelse(pvalue < 0.05,
                               ifelse(pvalue < 0.01, "**", "*"),"."), "")
  results.confound_dep <- cbind(parhat, se, pvalue, EC1)
  colnames(results.confound_dep) <- c("Estimate", "St.Dev.", "p", "CI.lb", "CI.ub")
  rownames(results.confound_dep) <- namescoef
  
  summary <- data.frame(round(results.confound_dep, 3))
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
  
  summary1 <- data.frame(round(results.naive, 3))
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
  
  summary2 <- data.frame(round(results.indep, 3))
  summary2$sign <- significant.indep
  summary2 <- summary2[,c(1:3, 6, 4:5)]
  summary2
  
  ## Create LaTeX tables of results
  xtab = xtable(summary)
  header= c("sample size",n,"Results 2-step_Estimation with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_2-step_Estimation_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary1)
  header= c("sample size",n,"Results naive model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_naive_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary2)
  header= c("sample size",n,"Results independence model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_independence_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  #
  # Plots of the estimated survival curves.
  #
  
  # parameter vector:
  # c(intercept, age, has_highschool_degree(1=yes),
  #   participated_in_study(no=0, otherwise=1),
  #   assigned_group(control=0, treatment = 1))
  
  # XandW variables: (intercept, age, has_highschool_degree(1=yes),
  #                   assigned_group(control=0, treatment = 1))
  dd <- c(1, 30, 0, 1) 
  
  # Z variable: participated_in_study(no=0, otherwise=1)
  Zobs <- 1
  
  #### Model with Yeo-Johnson transformation ####
  
  S <- NULL
  theta <- parhat[length(parhat)-1]
  s1 <- parhat[totparl + 1]
  gamma <- parhatG[(length(parhat)+1):length(parhatG)]
  
  V.obs <- (1-Zobs)*((1+exp(dd%*%gamma))*log(1+exp(dd%*%gamma))-(dd%*%gamma)*exp(dd%*%gamma))-Zobs*((1+exp(-(dd%*%gamma)))*log(1+exp(-(dd%*%gamma)))+(dd%*%gamma)*exp(-(dd%*%gamma))) #control function
  dd.2 <- c(dd[-length(dd)],Zobs,V.obs) #X,Z,V
  
  Time <- sort(exp(Y))
  
  for (i in 1:length(Time)) {
    sd = (YJtrans(log(Time[i]), theta) - t(parhat[1:length(dd.2)]) %*% dd.2)/s1
    S[i] = 1 - pnorm(sd)
  }
  
  #### Model without Yeo-Johnson transformation ####
  
  LikI_noTransform = function(par,Y,Delta,M) { # From Crommen, VK (2022)
    M=as.matrix(M)
    k = ncol(M)
    l = 2*k
    v = k+1
    beta = as.matrix(par[1:k])
    eta = as.matrix(par[v:l])
    sigma1 = par[l+1]
    sigma2 = par[l+2]
    
    z1 = (Y-(M%*%beta))/sigma1
    z2 = (Y-(M%*%eta))/sigma2
    
    tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta))
    p1 = pmax(tot,1e-100)
    Logn = sum(log(p1)); 
    return(-Logn)
  }
  
  LikF_noTransform = function(par,Y,Delta,M){ # From Crommen, VK (2022)
    M=as.matrix(M)
    k = ncol(M)
    l = 2*k
    v = k+1
    beta = as.matrix(par[1:k])
    eta = as.matrix(par[v:l])
    sigma1 = par[l+1]
    sigma2 = par[l+2]
    rho = par[l+3]
    
    z1 = (Y-(M%*%beta))/sigma1
    z2 = ((1-(rho*sigma2/sigma1))*Y-(M%*%eta-rho*(sigma2/sigma1)*(M%*%beta)))/(sigma2*((1-rho^2)^0.5))
    z3 = (Y-(M%*%eta))/sigma2
    z4 = ((1-(rho*sigma1/sigma2))*Y-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5)
    tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)
    p1 = pmax(tot,1e-100)   
    Logn = sum(log(p1)); 
    return(-Logn)
  }
  
  # Create matrix of X, Z and V.
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl),1,1) # Starting values
  
  # Independent model for starting values:
  parhat1 <- nloptr(x0=c(init),eval_f=LikI_noTransform,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5),ub=c(rep(Inf,totparl),Inf,Inf),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  initd = c(parhat1,0)
  parhat_noTransform <- nloptr(x0=initd,eval_f=LikF_noTransform,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99),ub=c(rep(Inf,totparl),Inf,Inf,0.99),
                               eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  
  
  S_noTransform <- NULL
  s1_noTransform <- parhat[totparl + 1]
  
  for (i in 1:length(Time)) {
    sd = (log(Time[i]) - t(parhat_noTransform[1:length(dd.2)]) %*% dd.2)/s1_noTransform
    S_noTransform[i] = 1 - pnorm(sd)
  }
  
  #### Plot ####
  
  print(paste("No Transformation: ",exp(qnorm(0.5)*s1_noTransform+t(parhat_noTransform[1:length(dd.2)]) %*% dd.2)))
  print(paste("Transformation: ",exp(IYJtrans(qnorm(0.5)*s1+t(parhat[1:length(dd.2)]) %*% dd.2,theta))))
  
  plot(Time, S, type = 'l', col = 1, xlab = "Time (t)", ylab = "S(t)",
       xlim = c(0, 365))
  lines(Time, S_noTransform, type = 'l', col = 1, lty=3)
  
  # Add legend
  legend("topright", c("Transformation model", "No transformation"),
         col = c(1, 1), lty = c(1,3))
}

########################### Data application Breast Cancer ##############################

DataApplicationBC <- function(data, init.value.theta_1, init.value.theta_2,
                              multiple.starting.points = FALSE) {
  
  # Remove this once the two "stop()" functions are resolved (somewhere around
  # code lines 4097 and 4246)
  #stop("Make sure all stops are removed from this code before fully running it")
  
  # Extract the parameters
  n = nrow(data)
  Y = data[, 1]
  Delta = data[, 2]
  Xi = data[, 3]
  intercept = data[,4]
  X = data[,(5:(parl + 1))]
  Z = data[,parl+2]
  W = data[,parl+3]
  XandW = as.matrix(cbind(data[,4],X,W))
  Xandintercept = cbind(intercept,X)
  
  # Estimate V
  gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma_BC,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  V <- (W-Z)*((1+exp(Xandintercept%*%gammaest))*log(1+exp(Xandintercept%*%gammaest))-(Xandintercept%*%gammaest)*exp(Xandintercept%*%gammaest))-Z*((1+exp(-(Xandintercept%*%gammaest)))*log(1+exp(-(Xandintercept%*%gammaest)))+(Xandintercept%*%gammaest)*exp(-(Xandintercept%*%gammaest)))
  
  # Create matrix of X, Z and V.
  M <- cbind(data[,4:(2+parl)],V)
  
  # Create matrix of X, Z and W.
  MnoV = data[,4:(3+parl)]
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  
  # Independent model for starting values sigmas and theta.
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  #
  # Our model: Taking into account Z is likely a confounded variable and that T 
  #            and C are dependent
  #
  
  # Create vector of initial values for the estimation of the parameters using the
  # parhat1. The final vector will be of the following form. 
  # [1:7]  : beta
  # [8:14] : eta
  # [15]   : sigma1
  # [16]   : sigma2
  # [17]   : rho
  # [18]   : theta_1
  # [19]   : theta_2
  
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
    
    initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
    initd[length(initd) - 2] <- 0
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
  }
  
  # Concatenate estimated coefficients vector for the models for transformed event
  # and censoring time with the vector of coefficients for the model
  # Z = \gamma * W + V. This results in the following vector.
  # [1:7]   : beta
  # [8:14]  : eta
  # [15]    : sigma1
  # [16]    : sigma2
  # [17]    : rho
  # [18]    : theta_1
  # [19]   : theta_2
  # [20:25] : gamma
  parhatG = c(parhat,as.vector(gammaest))
  
  Hgamma = hessian(LikFG_BC,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=10, v=2, show.details=FALSE)) 
  
  # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
  # (i.e. H_delta).
  H = Hgamma[1:length(initd),1:length(initd)]
  HI = ginv(H)
  
  Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
  
  prodvec = XandW[,1]
  
  for (i in 1:parlgamma) {
    for (j in 2:parlgamma) {
      if (i<=j){
        prodvec<-cbind(prodvec,diag(Xandintercept[,i]%*%t(Xandintercept[,j])))
      }
    }
  }
  
  secder=t(-dlogis(Xandintercept%*%gammaest))%*%prodvec
  
  WM = secder[1:parlgamma]
  WM <- rbind(WM, secder[c(2,3)])
  WMI = ginv(WM)
  
  diffvec = Z-plogis(Xandintercept%*%gammaest)*W
  
  mi = c()
  
  for(i in 1:n){
    newrow<-diffvec[i,]%*%Xandintercept[i,]
    mi = rbind(mi,newrow)
  }
  
  psii = -WMI%*%t(mi)
  
  gi = c()
  
  # For debugging
  #stop("Check whether this should be M[i,] or t(M[i,]) below")
  # In de LikF functie worden de kolommen van M[i,] geteld om de waarde van k
  # te bepalen. Als we dus bv. verkeerdelijk t(M[i,]) meegeven ipv M[i,], dan
  # is k = 1 in plaats van k = (iets anders). Bijgevolg worden dan de verkeerde
  # waarden gebruikt voor de parameters in die functie.
  
  # In de andere data applicatie staat er t(M[i,]), en is dit ook correct. Echter,
  # toen ik de code runde met mijn dummy data set ging het wel mis...
  
  for (i in 1:n) {
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
  
  
  #
  # Naive model: assuming Z is an unconfounded variable, but including dependence
  #              between T and C.
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
  
  #
  # Independence model: Model taking into account that Z is likely a confounded
  #                     variable but assuming independence between T and C.
  #
  
  # We construct the vector with
  # [1:7]   : beta
  # [8:14]  : eta
  # [15]    : sigma1
  # [16]    : sigma2
  # [17]    : theta_1
  # [18]    : theta_2
  # [19:24] : gamma
  
  parhatGI = c(parhat1,as.vector(gammaest))
  
  HgammaI = hessian(LikIGamma_BC,parhatGI,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  
  HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
  HIInd = ginv(HInd)
  
  VargammaI = HgammaI[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
  
  giI = c()
  
  # For debugging
  #stop("Check whether this should be M[i,] or t(M[i,])")
  # Zelfde opmerking als voordien.
  
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
  
  
  # Results of model assuming confounding and dependence between T and C.
  pvalue <- 2*pmin((1-pnorm(parhat/se)),pnorm(parhat/se))
  significant <- ifelse(pvalue < 0.10,
                        ifelse(pvalue < 0.05,
                               ifelse(pvalue < 0.01, "**", "*"),"."), "")
  results.confound_dep <- cbind(parhat, se, pvalue, EC1)
  colnames(results.confound_dep) <- c("Estimate", "St.Dev.", "p", "CI.lb", "CI.ub")
  rownames(results.confound_dep) <- namescoef
  
  summary <- data.frame(round(results.confound_dep, 3))
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
  
  summary1 <- data.frame(round(results.naive, 3))
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
  
  summary2 <- data.frame(round(results.indep, 3))
  summary2$sign <- significant.indep
  summary2 <- summary2[,c(1:3, 6, 4:5)]
  summary2
  
  ## Create LaTeX tables of results
  xtab = xtable(summary)
  header= c("sample size",n,"Results 2-step_Estimation with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_2-step_Estimation_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary1)
  header= c("sample size",n,"Results naive model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_naive_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  xtab = xtable(summary2)
  header= c("sample size",n,"Results independence model with YT-transformation")
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  # print.xtable(xtab,file=paste0("Results_independence_YT",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  
  
  #
  # Plots of the estimated survival curves.
  #
  
  # parameter vector:
  # c(intercept, age, participated_in_study(no=0, otherwise=1),
  #   invited_to_study(no=0, yes = 1))
  
  # XandW variables: (intercept, age, invited_to_study(no=0, yes = 1))
  dd <- c(1, -1.63, 1)
  
  # Z variable: participated_in_study(no=0, otherwise=1)
  Zobs <- 1
  
  #### Model with Yeo-Johnson transformation ####
  
  S <- NULL
  theta <- parhat[length(parhat)-1]
  s1 <- parhat[totparl + 1]
  gamma <- parhatG[(length(parhat)+1):length(parhatG)]
  
  V.obs <- (dd[3]-Zobs)*((1+exp(dd[-3]%*%gamma))*log(1+exp(dd[-3]%*%gamma))-(dd[-3]%*%gamma)*exp(dd[-3]%*%gamma))-Zobs*((1+exp(-(dd[-3]%*%gamma)))*log(1+exp(-(dd[-3]%*%gamma)))+(dd[-3]%*%gamma)*exp(-(dd[-3]%*%gamma))) #control function
  dd.2 <- c(dd[-length(dd)],Zobs,V.obs) #X,Z,V
  
  Time <- sort(exp(Y))
  
  for (i in 1:length(Time)) {
    sd = (YJtrans(log(Time[i]), theta) - t(parhat[1:length(dd.2)]) %*% dd.2)/s1
    S[i] = 1 - pnorm(sd)
  }
  
  #### Model without Yeo-Johnson transformation ####
  
  LikI_noTransform = function(par,Y,Delta,M) { # From Crommen, VK (2022)
    M=as.matrix(M)
    k = ncol(M)
    l = 2*k
    v = k+1
    beta = as.matrix(par[1:k])
    eta = as.matrix(par[v:l])
    sigma1 = par[l+1]
    sigma2 = par[l+2]
    
    z1 = (Y-(M%*%beta))/sigma1
    z2 = (Y-(M%*%eta))/sigma2
    
    tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta))
    p1 = pmax(tot,1e-100)
    Logn = sum(log(p1)); 
    return(-Logn)
  }
  
  LikF_noTransform = function(par,Y,Delta,M){ # From Crommen, VK (2022)
    M=as.matrix(M)
    k = ncol(M)
    l = 2*k
    v = k+1
    beta = as.matrix(par[1:k])
    eta = as.matrix(par[v:l])
    sigma1 = par[l+1]
    sigma2 = par[l+2]
    rho = par[l+3]
    
    z1 = (Y-(M%*%beta))/sigma1
    z2 = ((1-(rho*sigma2/sigma1))*Y-(M%*%eta-rho*(sigma2/sigma1)*(M%*%beta)))/(sigma2*((1-rho^2)^0.5))
    z3 = (Y-(M%*%eta))/sigma2
    z4 = ((1-(rho*sigma1/sigma2))*Y-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5)
    tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)
    p1 = pmax(tot,1e-100)   
    Logn = sum(log(p1)); 
    return(-Logn)
  }
  
  # Create matrix of X, Z and V.
  M = cbind(data[,4:(2+parl)],V)
  
  init = c(rep(0,totparl),1,1) # Starting values
  
  # Independent model for starting values:
  parhat1 <- nloptr(x0=c(init),eval_f=LikI_noTransform,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5),ub=c(rep(Inf,totparl),Inf,Inf),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  initd = c(parhat1,0)
  parhat_noTransform <- nloptr(x0=initd,eval_f=LikF_noTransform,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99),ub=c(rep(Inf,totparl),Inf,Inf,0.99),
                               eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  
  
  
  S_noTransform <- NULL
  s1_noTransform <- parhat[totparl + 1]
  
  for (i in 1:length(Time)) {
    sd = (log(Time[i]) - t(parhat_noTransform[1:length(dd.2)]) %*% dd.2)/s1_noTransform
    S_noTransform[i] = 1 - pnorm(sd)
  }
  
  #### Plot ####
  
  print(paste("No Transformation: ",exp(qnorm(0.5)*s1_noTransform+t(parhat_noTransform[1:length(dd.2)]) %*% dd.2)))
  print(paste("Transformation: ",exp(IYJtrans(qnorm(0.5)*s1+t(parhat[1:length(dd.2)]) %*% dd.2,theta))))
  
  plot(Time, S, type = 'l', col = 1, xlab = "Time (t)", ylab = "S(t)",
       xlim = c(0, max(Time)), ylim=c(0.96,1))
  lines(Time, S_noTransform, type = 'l', col = 1, lty=3)
  
  # Add legend
  legend("topright", c("Transformation model", "No transformation"),
         col = c(1, 1), lty = c(1,3))
}

########################## Data application chess ##############################

DataApplicationChess <- function(data, init.value.theta_1, init.value.theta_2) {
  
  
  n = nrow(data)
  Y = data[,1]
  Delta = data[,2]
  Xi = data[,3]
  intercept = data[,4]
  X = data[,5:(parl + 1)]
  Z = data[,parl+2]
  W = data[,parl+3]
  XandW = cbind(intercept,X,W)

  
  # Estimate V
  gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  ##############################################################################
  # This works better
  gammaest2 <- summary(glm(as.factor(Z) ~ -1 + XandW, family = "binomial"))$coefficients
  
  LikGamma2(gammaest, Z, XandW)
  LikGamma2(gammaest2[,1], Z, XandW)
  
  gammaest <- gammaest2[,1]
  ##############################################################################
  
  V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))

  # Create matrix of X, Z and V.
  M <- cbind(data[,4:(2+parl)],V)
  
  # Create matrix of X, Z and W.
  MnoV = data[,4:(3+parl)]
  
  init = c(rep(0,totparl), 1, 1, init.value.theta_1, init.value.theta_2)
  
  # Independent model for starting values sigmas and theta.
  parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5, 0,0),ub=c(rep(Inf,totparl),Inf,Inf, 2,2),
                   eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  #
  # Our model: Taking into account Z is likely a confounded variable and that T 
  #            and C are dependent
  #
  
  # Create vector of initial values for the estimation of the parameters using the
  # parhat1. The final vector will be of the following form. 
  # [1:7]  : beta
  # [8:14] : eta
  # [15]   : sigma1
  # [16]   : sigma2
  # [17]   : rho
  # [18]   : theta_1
  # [19]   : theta_2
  
  initd <-  c(parhat1[-length(parhat1)],parhat1[length(parhat1)-1],parhat1[length(parhat1)])
  initd[length(initd) - 2] <- 0
  
  parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,Xi=Xi,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-0.99,0,0),ub=c(rep(Inf,totparl),Inf,Inf,0.99,2,2),
                  eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
  
  # Concatenate estimated coefficients vector for the models for transformed event
  # and censoring time with the vector of coefficients for the model
  # Z = \gamma * W + V. This results in the following vector.
  # [1:7]   : beta
  # [8:14]  : eta
  # [15]    : sigma1
  # [16]    : sigma2
  # [17]    : rho
  # [18]    : theta_1
  # [19]   : theta_2
  # [20:25] : gamma
  parhatG = c(parhat,as.vector(gammaest))
  
  Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,Xi=Xi,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
  
  # Select part of variance matrix pertaining to beta, eta, var1, var2, rho and theta
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
  WM<-rbind(WM,secder[c(2,6,7,8,9)])
  WM<-rbind(WM,secder[c(3,7,10,11,12)])
  WM<-rbind(WM,secder[c(4,8,11,13,14)])
  WM<-rbind(WM,secder[c(5,9,12,14,15)])
  

  WMI = ginv(WM)
  
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
  
  
  #
  # Naive model: assuming Z is an unconfounded variable, but including dependence
  #              between T and C.
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
  
  #
  # Independence model: Model taking into account that Z is likely a confounded
  #                     variable but assuming independence between T and C.
  #
  
  # We construct the vector with
  # [1:7]   : beta
  # [8:14]  : eta
  # [15]    : sigma1
  # [16]    : sigma2
  # [17]    : theta_1
  # [18]    : theta_2
  # [19:24] : gamma
  
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
  
  
  # Results of model assuming confounding and dependence between T and C.
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
  
  return(c(parhat,gammaest))
}


