#### Example of an MCMC algorithm for fitting a location mixture of zif poisson components
#### The algorithm is tested using simulated data

## Suppose we have a mixture of 2 normal distribution.
##  The pdf of zif poisson is corresponding to:
##      f(x) = w * I(·) + (1-w) * Pois(lambda)
##  and we have 120 observations.

## Clear the environment and load required libraries
rm(list=ls())
#install.packages("MCMCpack")
library(MCMCpack)
set.seed(81196)  # So that results are reproducible


######## step1: load data #########
size=read.csv("nestsize.csv",stringsAsFactors = F)
X = size$X4
View(X)

## Plot the Empirical distribution
table(X)
par(mfrow=c(1,1))
barplot(table(X))
length(X[X ==0]) /length(X)

######## step2: Initialize the parameters for gibbs sampling #########
KK    =  2   # num of component
w     = length(X[X ==0]) /length(X) # weights of zero coef          
lambda    = 3   # guess of poisson expectation
n=length(X)

## Plot the initial guess for the density
xx = seq(0,11)
yy = w*c(1,rep(0,length(xx)-1)) + (1-w)*dpois(xx,lambda = lambda)
barplot(yy, ylim=c(0, max(yy)), xlab="x", ylab="Initial Mass")


## The actual MCMC algorithm starts here
# Priors
aa  = rep(1,KK)  # Uniform prior on w (params of beta dist prior on w)
lambda_0 = 1          # exp prior for lambda 

# Number of iterations of the sampler
rrr   = 5000
burn  = 1000


# Storing the samples
## used to store the indicators(where xi comes from)
## an array of 2-dim (same as dataframe),row-ith contains the assignment of ith iteration.
## col-jth contains the jth rank of sample.
cc.out    = array(0, dim=c(rrr, n))  
w.out     = rep(0, rrr)                 
lambda.out    = rep(0, rrr)   
logpost   = rep(0, rrr)  ## logarithm of posterior likelihood, used to monitor the convergence of algorithm

######## step3: MCMC iterations ##########
for(s in 1:rrr){
  # Sample the indicators c_i for each sample x_i(with a for loop)
  cc = rep(0,n)
  for(i in 1:n){
    v = rep(0,KK)
    v[1] = w * ifelse(X[i]>0,0,1)
    #  v[2] = (1-w) * dpois(X[i],lambda = lambda) , both OK, since we normalized after (v = v/sum(v))
    v[2] = (1-w) * ifelse(X[i]>0,1,dpois(X[i],lambda = lambda))
    v = v/sum(v)
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v) ## category
  }
  
  # Sample the weights
  w = rbeta(1, 1 + sum(cc==1), 1 + sum(cc==2)) ## uniform prior
  
  # Sample the lambda
  lambda = rgamma(1, 1 + sum(X[cc==2]), lambda_0 + sum(cc==2))
  
  # Store samples
  cc.out[s,]   = cc
  w.out[s]     = w
  lambda.out[s]   = lambda
  
  # X
  for(i in 1:n){
    if(cc[i]==1){
      logpost[s] = logpost[s] + log(w)
    }else{
      logpost[s] = logpost[s] + log(1-w) + dpois(X[i],lambda = lambda, log=TRUE)
    }
  }
  # w
  logpost[s] = logpost[s] + dbeta(w, aa[1], aa[2],log=T) ## log=T?
  # lambda
  logpost[s] = logpost[s] + dgamma(lambda, 1, lambda_0)
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}


## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")


## posterior estimate of lambda, weights
round(mean(lambda.out[burn:rrr]),2) # 3.05
round(mean(w.out[burn:rrr]),2) # 0.4
