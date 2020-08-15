#### Semisupervised classification for the banknote dataset
rm(list=ls())
library(mvtnorm)
library(MCMCpack)

## Load data
load("banknoteclassification.Rdata")
x = rbind(banknote.training,banknote.test)

## Generate data from a mixture with 3 components
KK      = length(unique(banknote.training.labels))
p       = dim(banknote.training)[2]
n       = dim(banknote.training)[1]
m       = dim(unique(banknote.test))[1]

## Initialize the parameters
w          = rep(1,KK)/KK  #Assign equal weight to each component to start with
mu         = rmvnorm(KK, apply(x,2,mean), var(x))   #RandomCluster centers randomly spread over the support of the data
Sigma      = array(0, dim=c(KK,p,p))  #Initial variances are assumed to be the same
Sigma[1,,] = var(x)/KK  
Sigma[2,,] = var(x)/KK
cc         = c(as.numeric(banknote.training.labels), sample(1:KK, m, replace=TRUE, prob=w))

# Priors
aa = rep(1, KK)
dd = apply(x,2,mean)
DD = 10*var(x)
nu = p+1
SS = var(x)/3

# Number of iteration of the sampler
rrr  = 11000
burn = 1000

# Storing the samples
cc.out    = array(0, dim=c(rrr, n+m))
w.out     = array(0, dim=c(rrr, KK))
mu.out    = array(0, dim=c(rrr, KK, p))
Sigma.out = array(0, dim=c(rrr, KK, p, p))
logpost   = rep(0, rrr)

for(s in 1:rrr){
  # Sample the indicators
  for(i in (n+1):(n+m)){
    v = rep(0,KK)
    for(k in 1:KK){
      v[k] = log(w[k]) + dmvnorm(x[i,], mu[k,], Sigma[k,,], log=TRUE)  #Compute the log of the weights
    }
    v = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }
  
  # Sample the weights
  w = as.vector(rdirichlet(1, aa + tabulate(cc)))
  
  # Sample the means
  DD.st = matrix(0, nrow=p, ncol=p)
  for(k in 1:KK){
    mk    = sum(cc==k)
    xsumk = apply(x[cc==k,], 2, sum)
    DD.st = solve(mk*solve(Sigma[k,,]) + solve(DD))
    dd.st = DD.st%*%(solve(Sigma[k,,])%*%xsumk + solve(DD)%*%dd)
    mu[k,] = as.vector(rmvnorm(1,dd.st,DD.st))
  }
  
  # Sample the variances
  xcensumk = array(0, dim=c(KK,p,p))
  for(i in 1:(n+m)){
    xcensumk[cc[i],,] = xcensumk[cc[i],,] + (x[i,] - mu[cc[i],])%*%t(x[i,] - mu[cc[i],])
  }
  for(k in 1:KK){
    Sigma[k,,] = riwish(nu + sum(cc==k), SS + xcensumk[k,,])
  }
  
  # Store samples
  cc.out[s,]      = cc
  w.out[s,]       = w
  mu.out[s,,]     = mu
  Sigma.out[s,,,] = Sigma
  for(i in 1:n){
    logpost[s] = logpost[s] + log(w[cc[i]]) + dmvnorm(x[i,], mu[cc[i],], Sigma[cc[i],,], log=TRUE)
  }
  logpost[s] = logpost[s] + ddirichlet(w, aa)
  for(k in 1:KK){
    logpost[s] = logpost[s] + dmvnorm(mu[k,], dd, DD)
    logpost[s] = logpost[s] + diwish(Sigma[k,,], nu, SS)
  }
  
  if(s/250==floor(s/250)){
    print(paste("s = ", s))
  }
}
### qda
banknote.qda = qda(x=banknote.training, grouping=banknote.training.labels)
predict(banknote.qda, banknote.test)$class

### clf
probgenuine = rep(NA, m)
for(i in 1:m){
  probgenuine[i] = sum(cc.out[-seq(1,burn),n+i]==2)/(rrr-burn)
}
probgenuine[probgenuine >0.5] = 2
probgenuine[probgenuine <0.5] = 1
sum(probgenuine != as.numeric(banknote.test.labels))


### post estimate
mu1.out.df = as.data.frame(mu.out[burn:rrr,1,])
mu2.out.df = as.data.frame(mu.out[burn:rrr,2,])
Sigma1.out.df = as.data.frame(Sigma.out[burn:rrr,1,,])
Sigma2.out.df = as.data.frame(Sigma.out[burn:rrr,2,,])

mu1.post = apply(mu1.out.df,2,mean)
mu2.post = apply(mu2.out.df,2,mean)
# mu1.post
#       V1        V2        V3        V4        V5        V6 
# 214.82289 130.29968 130.19279  10.53259  11.13298 139.44921 
# mu2.post
#       V1         V2         V3         V4         V5         V6 
# 214.968776 129.943345 129.720413   8.305658  10.168912 141.516382 
mu1.post
mu2.post

Sigma1.post = matrix(apply(Sigma1.out.df,2,mean),p,p)
Sigma2.post = matrix(apply(Sigma2.out.df,2,mean),p,p)

Sigma1.post
Sigma2.post
