#### Example of an MCMC algorithm for fitting a location mixture of exponential and log-gaussian components
## Suppose we have a mixture of exponential and log-gaussian. 
##  The pdf of is corresponding to:
##      f(x) = w * exp(lambda) + (1-w) * lognormal(mu,tau)
##  and we have 120 observations.




## Clear the environment and load required libraries
rm(list=ls())
#install.packages("MCMCpack")
library(MCMCpack)
#set.seed(81196)  # So that results are reproducible


######## step1: load data #########
fuse=read.csv("fuses.csv",stringsAsFactors = F,header = F)
#View(fuse)
x = fuse$V1

## Plot the Empirical distribution
par(mfrow=c(1,1))
hist(x,breaks=50)


######## step2: Initialize the parameters for gibbs sampling #########
KK=2            # number of cluster
w = 0.1        # weight of expo
lambda    = 30   #param of Exponential-distribution.
mu = mean(log(x))    
tau = sd(log(x)) 
n=length(x)
cc = sample(1:2, n, replace=TRUE, prob=c(w,1-w))

## Plot the initial guess for the density
xx = seq(0,11,length=200)
yy = w*dexp(xx,rate = lambda) + (1-w)*dlnorm(xx,meanlog = mu,sdlog = sqrt(tau2))
plot(xx, yy, type="l", ylim=c(0, max(yy)), xlab="x", ylab="Initial density")


## The actual MCMC algorithm starts here
# Priors
aa  = rep(1,KK)  # Uniform prior on w (params of beta dist prior on w)
alpha = 1        # gamma prior for expo lambda  
beta =  1        # gamma prior for expo lambda  
mu_0 = 0          # normal prior for mu
tau_0 = 1        # normal prior for mu
d = 2            # inv-gamma prior for tau2
q = 1            # inv-gamma prior for tau2
# Number of iterations of the sampler
rrr   = 10000
burn  = 3000


# Storing the samples
## used to store the indicators(where xi comes from)
## an array of 2-dim (same as dataframe),row-ith contains the assignment of ith iteration.
## col-jth contains the jth rank of sample.
cc.out    = array(0, dim=c(rrr, n))  
w.out     = rep(0, rrr)                 
mu.out    = rep(0, rrr)   
tau.out  = rep(0, rrr)   
lambda.out  = rep(0, rrr)  
logpost   = rep(0, rrr)  ## logarithm of posterior likelihood, used to monitor the convergence of algorithm


######## step3: MCMC iterations ##########
#View(cc.out)
for(s in 1:rrr){
  # Sample the indicators c_i for each sample x_i(with a for loop)
  v = rep(0,KK)
  for(i in 1:n){
    v[1]  = log(w) + dexp(x[i], lambda, log=TRUE)
    v[2] = log(1-w) + dlnorm(x[i],mu,tau, log = TRUE)
    #Go from logs to actual weights in a numerically stable manner
    v     = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:2, 1, replace=TRUE, prob=v)## category  
  }
  
  # Sample the weights
  w = rbeta(1, aa[1] + sum(cc==1), aa[2] + sum(cc==2))
  
  # Sample the lambda
  lambda = rgamma(1, alpha + sum(cc==1), beta + sum(x[cc==1]))
  
  # Sample the mu
  n2    = sum(cc==2)
  mu.post = ( sum(log(x[cc==2])) / tau^2
             ) / ( (n2/tau^2) +1/tau_0^2 )
  sigma2.post = 1/(n2/tau^2 + 1/tau_0^2)
  mu   = rnorm(1, mu.post, sqrt(sigma2.post))
  
  # Sample the tau
  tau = sqrt(MCMCpack::rinvgamma(1,d + sum(cc==2), q + sum((log(x[cc==2]) -mu)^2)))

  # Store samples
  cc.out[s,]   = cc
  w.out[s]     = w
  lambda.out[s]   = lambda
  mu.out[s] = mu 
  tau.out[s] = tau
  # x
  for(i in 1:n){
    if(cc[i]==1){
      logpost[s] = logpost[s] + log(w) + dexp(x[i],rate = lambda,log=T)
    }else{
      logpost[s] = logpost[s] + log(1-w) + dlnorm(x[i],meanlog= mu,sdlog =tau, log=TRUE)
    }
  }
  # w
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2])
  # mu
  logpost[s] = logpost[s] + dnorm(mu, mu_0, tau_0) 
  # tau
  logpost[s] = logpost[s] + log(dinvgamma(tau, d, q)) 
  # lambda
  logpost[s] = logpost[s] + dgamma(lambda, 1, 1)
  
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}


## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")

##  "w=0.1 lamdba=2.29 mu=0.79 tau=0.38"
## posterior estimate of lambda, weights
round(mean(w.out[burn:rrr]),2) # 0.1
round(mean(mu.out[burn:rrr]),2) # 0.79
round(mean(tau.out[burn:rrr]),2) # 0.38
round(mean(lambda.out[burn:rrr]),2) # 2.29

