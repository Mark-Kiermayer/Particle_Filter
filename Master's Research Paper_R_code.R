#######################################################################
############## Master's Research Paper  ###############################
##############      Mark Kiermayer      ###############################
##############      Implementation      ###############################
#######################################################################


### Define parameters (see De Rossi, i.e. estimates by de Jong and Santa-Clara (1999). "The dynamics of the forward interest rate curve: A formulation with state variables" )
T = 150
kappa = 0.1862
theta = 0.0654
sigma = 0.0481
lambda = -32.03*sigma^2
N = 5
h = 0.01^2 # stdev of 10 basis points for measurement errors <- h = 0.001^2 results in normalizing constants K_i=0 -> weights = 0 -> statistic for resampling stat = NaN ..
dt = 1 # weekly observations, estimates for parameters by de Jong and Santa-Clara are based on weekly observations
N_sample = 100 # number of sample paths
# Remark: We create 100 paths to be able to evaluate the algorithms accuracy and variance of estimates

# Allow for run-time analysis of likelihood calculation
library(tictoc)

### Simulate short rates
alpha = matrix(NA,N_sample,T)
q = 2*kappa*theta/sigma^2-1
c = 2*kappa/(sigma^2*(1-exp(-kappa*dt)))

###########################################################################
#
# Simulation of Data
#
###########################################################################

# set seet for simmulations
set.seed(42)

# Initial values for CIR model, i.e. drawn from Gamma dist.
alpha[,1]= rgamma(N_sample, rate = c*(1-exp(-kappa*dt)),shape = q+1)
# Next we draw from conditional transition density, i.e. noncentral chi square
# Theory: non-centreal chi square is mixture of central chi squares with df proportional to poisson distr.
# Method: draw df from poisson dist. and then simulate chi square
# Method outlined by Duan and Simonate (1999)
for(i in 2:T){
  df = 2*q+2+2*rpois(N_sample, c*alpha[,i-1]*exp(-kappa*dt))
  g = rchisq(n=N_sample,df)
  alpha[,i] = g/(2*c)
}

# figure with selected paths
par(mfrow=c(1,1))
plot(alpha[1,], type='l', ylim = c(min(alpha),max(alpha)),xlab = 'week', ylab = "short rate", main = "CIR model - selected paths")
for(i in 1:10){
  lines(alpha[i,], col = i)
}
abline(h=theta, lty = 2, lwd = 2)
legend('topleft', legend = 'long-term mean', lty =2)

# Simulate bonds, i.e. yields, for maturities of 0.25, 1, 3 and 5 for each path of simulated short rates alpha
y_025 = matrix(NA,N_sample,T)
y_1 = matrix(NA,N_sample,T)
y_3 = matrix(NA,N_sample,T)
y_5 = matrix(NA,N_sample,T)
y_10 = matrix(NA,N_sample,T)
gamma = sqrt((kappa+lambda)^2+2*sigma^2)

# Functions for affine yield curve
Yield_A <- function(tau,kappa, lambda, gamma, theta){
  A = (2*kappa*theta/sigma^2)/tau*log((2*gamma*exp(tau*(kappa+lambda+gamma)/2))/(2*gamma+(kappa+lambda+gamma)*(exp(gamma*tau)-1)))
  return(A)
}
Yield_B <- function(tau, kappa, lambda, gamma){
  B = 2/tau*(exp(gamma*tau)-1)/(2*gamma+(kappa+lambda+gamma)*(exp(gamma*tau)-1))
  return(B)
}

# Note: we have Gaussian iid error term, i.e. we can simulate them individually
for(i in 1:T){
  y_025[,i] = -Yield_A(0.25,kappa, lambda, gamma, theta)+Yield_B(0.25,kappa, lambda, gamma)*alpha[,i]+rnorm(N_sample,0,sqrt(h))
  y_1[,i] = -Yield_A(1,kappa, lambda, gamma, theta)+Yield_B(1,kappa, lambda, gamma)*alpha[,i]+rnorm(N_sample,0,sqrt(h))
  y_3[,i] = -Yield_A(3,kappa, lambda, gamma, theta)+Yield_B(3,kappa, lambda, gamma)*alpha[,i]+rnorm(N_sample,0,sqrt(h))
  y_5[,i] = -Yield_A(5,kappa, lambda, gamma, theta)+Yield_B(5,kappa, lambda, gamma)*alpha[,i]+rnorm(N_sample,0,sqrt(h))
  y_10[,i] = -Yield_A(10,kappa, lambda, gamma, theta)+Yield_B(10,kappa, lambda, gamma)*alpha[,i]+rnorm(N_sample,0,sqrt(h))
}

############################################################################
#
# Main function - Log-likelihood estimation using particle filter method
#
############################################################################

library(VGAM)
# Calculate K_i explicitely, check Wolfram Alpha
K_int <- function(a,b,c, N,h){
  const = 1/(2*pi*h)^(N/2)
  int = sqrt(pi)*exp(b^2/(4*a)-c)*erfc(b/(2*sqrt(a)))/(2*sqrt(a))
  return(const*int)
}

library(truncnorm)
# function for particle filter estimation of likelihood
LH_estimate <- function(y_025,y_1,y_3,y_5,y_10,lambda_PF,kappa_PF, theta_PF, sigma_PF, N_part,h,dt, resample ){
  
  # derive resulting parameters
  gamma_PF = sqrt((kappa_PF+lambda_PF)^2+2*sigma_PF^2)
  c_PF = 2*kappa_PF/(sigma_PF^2*(1-exp(-kappa_PF*dt)))
  q_PF = 2*kappa_PF*theta_PF/sigma^2-1
  T = length(y_025)
  
  # store particles
  part = matrix(0,N_part,T)
  # estimates K_i
  K = matrix(0,T,1)
  # weights
  w = matrix(0,N_part,T)
  # count resampling steps
  count_resample = 0
  # parameters for importance density
  b_t = matrix(0,1,T)
  c_t = matrix(0,1,T)
  # additional parameters for importance sampling, i.e. mean and variance of Gaussian importance density
  a = 1/(2*h)*(Yield_B(0.25,kappa_PF, lambda_PF, gamma_PF)^2+Yield_B(1,kappa_PF, lambda_PF, gamma_PF)^2+Yield_B(3,kappa_PF, lambda_PF, gamma_PF)^2+Yield_B(5,kappa_PF, lambda_PF, gamma_PF)^2+Yield_B(10,kappa_PF, lambda_PF, gamma_PF)^2)
  b_t[1] = -1/h*((y_025[1]-Yield_A(0.25,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(0.25,kappa_PF, lambda_PF, gamma_PF)+(y_1[1]-Yield_A(1,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(1,kappa_PF, lambda_PF, gamma_PF)+(y_3[1]-Yield_A(3,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(3,kappa_PF, lambda_PF, gamma_PF)+(y_5[1]-Yield_A(5,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(5,kappa_PF, lambda_PF, gamma_PF)+(y_10[1]-Yield_A(10,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(10,kappa_PF, lambda_PF, gamma_PF))
  c_t[1] = 1/(2*h)*((y_025[1]-Yield_A(0.25,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_1[1]-Yield_A(1,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_3[1]-Yield_A(3,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_5[1]-Yield_A(5,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_10[1]-Yield_A(10,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2)
  
  # simulate particles for first step (including truncation)
  IS_mean = -b_t[1]/(2*a)
  IS_var = 1/(2*a)
  part[,1]= rtruncnorm(N_part,0,Inf,rep(IS_mean,N_part),rep(sqrt(IS_var),N_part))
  # calculate K_1
  K[1] = K_int(a,b_t[1],c_t[1],N,h)
  # calculate weights
  w[,1] = K[1]*dgamma(part[,1],rate = c_PF*(1-exp(-kappa_PF*dt)), shape = q_PF+1)
  
  for(t in 2:T){
    
    # new parameter for IS
    b_t[t] = -1/h*((y_025[t]-Yield_A(0.25,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(0.25,kappa_PF, lambda_PF, gamma_PF)+(y_1[t]-Yield_A(1,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(1,kappa_PF, lambda_PF, gamma_PF)+(y_3[t]-Yield_A(3,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(3,kappa_PF, lambda_PF, gamma_PF)+(y_5[t]-Yield_A(5,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(5,kappa_PF, lambda_PF, gamma_PF)+(y_10[t]-Yield_A(10,kappa_PF, lambda_PF, gamma_PF, theta_PF))*Yield_B(10,kappa_PF, lambda_PF, gamma_PF))
    c_t[t] = 1/(2*h)*((y_025[t]-Yield_A(0.25,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_1[t]-Yield_A(1,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_3[t]-Yield_A(3,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_5[t]-Yield_A(5,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2+(y_10[t]-Yield_A(10,kappa_PF, lambda_PF, gamma_PF, theta_PF))^2)
    
    # simulate particles
    IS_mean = -b_t[t]/(2*a)
    #IS_var = 1/(2*a)
    part[,t]= rtruncnorm(N_part,0,Inf,rep(IS_mean,N_part),rep(sqrt(IS_var),N_part))
    # calculate K_j
    K[t] = K_int(a,b_t[t],c_t[t],N,h)
    # calculate weights
    trans = 2*c_PF*dchisq(2*c_PF*part[,t],df=2*q_PF+2,ncp=2*c_PF*exp(-kappa_PF*dt)*part[,t-1])
    
    # begin error control sequence
    if(any(trans==0)){
      print('error trans')
      print(t)
      return(part)    
    }
    # end error control sequence
    
    # iteratively determine new weights
    w[,t] = K[t]*trans*w[,t-1]/sum(w[,t-1])
    # Calculate statistic to decide whether to resample or not
    stat = sum(w[,t])^2/sum(w[,t]^2)
    
    # begin error control
    if(is.na(stat)==TRUE){
      print('error stat')
      return(list(weight=w[,], sim = trans,const =K, part = part))
    }
    # end error control
    
    # resampling of weights to prevent degeneracy
    if(stat<N_part/2&&resample==TRUE){ # resample is an input switch which allows to run the algorithm without resampling
      
      # keep track of number of resampling steps
      count_resample = count_resample +1
      
      cdf = cumsum(w[,t]/sum(w[,t])) # drop cdf_0=0 since not required and for the sake of indices
      cdf_hat = rep(1/N_part,N_part)
      cdf_hat[1] = runif(1,0,1/N_part)
      cdf_hat = cumsum(cdf_hat) # new sample cdf
      
      k = 1
      l = 1
      while(k<=N_part && l<=N_part){
        if(cdf[k]>=cdf_hat[l]){
          part[l,1:t] = part[k,1:t]
          l = l+1
        }else{
          k = k+1
        }
      }
      # adapt weights, i.e. equally weighted particles after resampling
      w[,t] = sum(w[,t])/N_part
    } # end of potential resampling algorithm  
  }
  LH_estimate = log(mean(w[,1]))+ sum(log(apply(w[,-1],2,sum)))
  return(list(LH = LH_estimate, weights = w, resample = count_resample))
}

############################################################################
#
# Let's use the simulated data and previously coded functions to perform a empirical analysis 
# In more detail I will show the importance of resampling, check the algorithm 

# Number of particles for MC
N_part = 100
# Parameters for particle filter, true parameters unknown (h assumed to be known)
lambda_PF = lambda
kappa_PF = kappa
theta_PF = theta
sigma_PF = sigma
gamma_PF = sqrt((kappa_PF+lambda_PF)^2+2*sigma_PF^2)
c_PF = 2*kappa_PF/(sigma_PF^2*(1-exp(-kappa_PF*dt)))
q_PF = 2*kappa_PF*theta_PF/sigma^2-1

#######################################
######## Effect of resampling #########
#######################################

# Run PF estimation for likelihood for single paths to observe necessity of resampling
output_no_RS = LH_estimate(y_025[1,],y_1[1,],y_3[1,],y_5[1,],y_10[1,],lambda_PF,kappa_PF, theta_PF, sigma_PF, N_part, h,dt, resample=FALSE )
weights_no_RS = matrix(unlist(output_no_RS$weights),N_part,T)

# Plot to highlight strong degeneracy if we DO NOT resample
par(mfrow=c(1,2))
plot(1:T,apply(weights_no_RS,2,max)/colSums(weights_no_RS[,1:T]),xlab = "time", ylab = "relative weight of dominating particle", type = 'l', main = 'Degeneracy - No Resampling', cex.main= 0.75)

# Now let's include the resampling algorithm
output_with_RS = LH_estimate(y_025[1,],y_1[1,],y_3[1,],y_5[1,],y_10[1,],lambda_PF,kappa_PF, theta_PF, sigma_PF, N_part, h,dt, resample=TRUE)
weights_with_RS = matrix(unlist(output_with_RS[2]),N_part,T)
count_resample = unlist(output_with_RS[3])
count_resample

# Plot to highlight diminishing degeneracy if we DO resample
plot(1:T,apply(weights_with_RS,2,max)/colSums(weights_with_RS[,1:T]),xlab = "time", ylab = "relative weight of dominating particle", type = 'l', main = 'Degeneracy - With Resampling', cex.main= 0.75)

############################################################
###### Analysis of the particle filter's stability #########
############################################################

# For now, let's use the true (in practice unknown) parameters.
# A parameter estimation will be the next step
## Average running time, avg. number of resampling steps and avg. std.dev. of log-likelihood for different number of particles

### 100 particles ###
output_100 = matrix(0,N_sample,1)
time_100 = matrix(0,N_sample,1)
RS_100 = matrix(0,N_sample,1)

# test procedure (using same path)
for(j in 1:N_sample){
  tic()
  output_100 = LH_estimate(y_025[1,],y_1[1,],y_3[1,],y_5[1,],y_10[1,],lambda_PF,kappa_PF, theta_PF, sigma_PF, N_part= 100,h,dt, resample=TRUE )
  x=toc(quiet = TRUE)
  LH_100[j] = unlist(output_100$LH)
  time_100[j]=x$toc-x$tic
  RS_100[j] = unlist(output_100$resample)
}

# Standard deviation of LH procedure
sqrt(var((LH_100)))
# average running time
mean(time_100)
#average resampling steps
mean(RS_100)

#### Now let's increase the number of simulated particles ###
### 200 Particles ###
output_200 = matrix(0,N_sample,1)
time_200 = matrix(0,N_sample,1)
RS_200 = matrix(0,N_sample,1)

# test procedure (using same path)
for(j in 1:N_sample){ 
  tic()
  output_200 = LH_estimate(y_025[1,],y_1[1,],y_3[1,],y_5[1,],y_10[1,],lambda_PF,kappa_PF, theta_PF, sigma_PF, N_part= 200,h,dt, resample=TRUE )
  x=toc(quiet = TRUE)
  LH_200[j] = unlist(output_200$LH)
  time_200[j]=x$toc-x$tic
  RS_200[j] = unlist(output_200$resample)
}

# Standard deviation of LH procedure
sqrt(var(LH_200))
# average running time
mean(time_200)
#average resampling steps
mean(RS_200)

### 500 Particles ###
output_500 = matrix(0,N_sample,1)
time_500 = matrix(0,N_sample,1)
RS_500 = matrix(0,N_sample,1)

# test procedure (using same path)
for(j in 1:N_sample){
  tic()
  output_500 = LH_estimate(y_025[1,],y_1[1,],y_3[1,],y_5[1,],y_10[1,],lambda_PF,kappa_PF, theta_PF, sigma_PF, N_part= 500,h,dt, resample=TRUE )
  x=toc(quiet = TRUE)
  LH_500[j] = unlist(output_500$LH)
  time_500[j]=x$toc-x$tic
  RS_500[j] = unlist(output_500$resample)
}

# Standard deviation of LH procedure
sqrt(var(LH_500))
# average running time
mean(time_500)
#average resampling steps
mean(RS_500)

###########################################################################
########### Parameter estimation using PF and LH approach #################
###########################################################################

# determine parameters by PF likelihood approach
# Let's apply a one-factor estimation to visualize the estimator
# Also, let's assume all true parameters except 'theta' to be known
# Let's use a grid search to look for a optimal parameter using our PF-LH approach

##### Estimate theta, all other parameters fixed
theta_PF_grid = sort(c(seq(0.005,0.2, 0.005), theta)) # grid - include true parameter
n_grid = length(theta_PF_grid)
LH_theta_est = matrix(0,N_sample,n_grid)
for(i in 1:N_sample){ # Use multiple simulations to allow for some variability of the estimation procedure and eventually use their mean as an estimate
  for(j in 1:n_grid){
    LH_theta_est[i,j] = unlist(LH_estimate(y_025[i,],y_1[i,],y_3[i,],y_5[i,],y_10[i,],lambda_PF,kappa_PF, theta_PF_grid[j], sigma_PF, N_part,h,dt, resample=TRUE )$LH)
  }
}
x = apply(LH_theta_est,1,which.max)
theta_est = theta_PF_grid[x]
theta_est
theta_PF_dist =hist(theta_est, plot = FALSE, breaks = n_grid)
theta_PF_dist$density = theta_PF_dist$counts/N_sample
par(mfrow = c(1,1))
plot(theta_PF_dist, freq = FALSE, main = 'distribution of estimate', xlim = c(0.00,0.2))
abline(v = theta, col = 'green')
abline(v = mean(theta_est), col = 'red')
legend('topright',legend = c("LH-estimate","estimated theta","true theta"), col = c("black","orange", "green"), lty =1)

## Plot Likelihood values for different values of theta ('N_sample' samples used)
par(mfrow=c(1,1))
plot(theta_PF_grid, LH_theta_est[1,], type='l', xlab = 'theta', ylab = 'Log-likelihood values', ylim = c(250,700), main ='Selected log-likelihood curves')
for(i in 2:(N_sample/5)){
  lines(theta_PF_grid,LH_theta_est[i,])#, col = i)
}
abline(v = theta, col = 'green')
legend('topright',legend = c("LH-estimate","true theta"), col = c("black", "green"), lty =1)