###### Setup data
x = faithful$eruptions
n = length(x)

### Get a "Bayesian" kernel density estimator based on the same location mixture of 6 normals
## Priors set up using an "empirical Bayes" approach
KK=2
aa  = rep(1,KK)  
eta = mean(x)    
tau = sqrt(var(x))
dd  = 2 
qq  = var(x)/KK

## Initialize the parameters
w     = 0.5
mu    = rnorm(KK, mean(x), sd(x))
sigma = rep(sd(x)/KK,KK)
cc    = sample(1:KK, n, replace=T, prob=c(w,1-w))

## Number of iterations of the sampler
rrr   = 1200
burn  = 200

## Storing the samples
cc.out    = array(0, dim=c(rrr, n))
w.out     = rep(0, rrr)
mu.out    = array(0, dim=c(rrr, KK))
sigma.out = array(0, dim=c(rrr, KK))
logpost   = rep(0, rrr) 



for(s in 1:rrr){
  # Sample the indicators
  cc = rep(0,n)
  for(i in 1:n){
    v = rep(0,KK)
    v[1] = log(w) + dnorm(x[i], mu[1], sigma[1], log=TRUE)  #Compute the log of the weights
    v[2] = log(1-w) + dnorm(x[i], mu[2], sigma[2], log=TRUE)  #Compute the log of the weights
    v = exp(v - max(v))/sum(exp(v - max(v)))
    #print(v)
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
    
  }
  
  # Sample the weights
  tmp = tabulate(cc, nbins=KK)
  w = rbeta(1,aa[1] + tmp[1],aa[2] + tmp[2])
  
  # Sample the means
  for(k in 1:KK){
    nk    = sum(cc==k)
    xsumk = sum(x[cc==k])
    tau2.hat = 1/(nk/sigma[k]^2 + 1/tau^2)
    mu.hat  = tau2.hat*(xsumk/sigma[k]^2 + eta/tau^2)
    mu[k]   = rnorm(1, mu.hat, sqrt(tau2.hat))
  }
    
  # Sample the variances
  for(k in 1:KK){
    nk    = sum(cc==k)
    dd.star = dd + nk/2
    qq.star = qq + sum((x[cc==k] - mu[k])^2)/2
    sigma[k] = sqrt(1/rgamma(1, dd.star, qq.star))
  }
  

  # Store samples
  cc.out[s,]   = cc
  w.out[s]    = w
  mu.out[s,]   = mu
  sigma.out[s,] = sigma
  for(i in 1:n){
  if(cc[i] == 1){
    logpost[s] = logpost[s] + log(w) + dnorm(x[i], mu[1], sigma[2], log=TRUE)
  } else{
   logpost[s] = logpost[s] + log(1-w) + dnorm(x[i], mu[2], sigma[2], log=TRUE)
    }
  }
  logpost[s] = logpost[s] + log(dbeta(w, aa[1],aa[2]))
  for(k in 1:KK){
    logpost[s] = logpost[s] + dnorm(mu[k], eta, tau, log=TRUE)
  }
  for(k in 1:KK){
    logpost[s] = logpost[s] + dgamma(1/sigma[k]^2, dd, qq, log=TRUE) - 4*log(sigma[k])
  }
  if(s/500==floor(s/500)){
    print(paste("s =",s))
  }
}
### posterior
round(mean(w.out[burn:rrr]),2) # 0.65
round(mean(mu.out[burn:rrr]),2) # 4.28

## Compute the samples of the density over a dense grid
xx  = seq(0,7,length=150)
density.mcmc = array(0, dim=c(rrr-burn,length(xx)))
for(s in 1:(rrr-burn)){
    density.mcmc[s,] = density.mcmc[s,] + 
    w.out[s+burn]*dnorm(xx,mu.out[s+burn,1],sigma.out[s+burn,1]) + 
    (1-w.out[s+burn])*dnorm(xx,mu.out[s+burn,2],sigma.out[s+burn,2])
}
density.mcmc.m = apply(density.mcmc , 2, mean)

length(xx)
length(density.mcmc.m)
density.mcmc.lq = apply(density.mcmc, 2, quantile, 0.025)
density.mcmc.uq = apply(density.mcmc, 2, quantile, 0.975)
plot(xx, density.mcmc.m, type="n",ylim=c(0,max(density.mcmc.uq)),xlab="", ylab="Density")
polygon(c(xx,rev(xx)), c(density.mcmc.lq, rev(density.mcmc.uq)), col="grey", border="grey")
lines(xx, density.mcmc.m, col=colscale[1], lwd=2)
points(x, rep(0,n))

