load("banknoteclassification.Rdata")
library(mvtnorm)
library(ellipse)
library(MCMCpack)

### inspect data ###
View(banknote.training) ## 30 rows, 6 features.
banknote.training.labels
View(banknote.test) ## 170 rows, 6 features.
banknote.test.labels


############ semi-supervised classification using MCMC ##################
set.seed(63252)    #Keep seed the same so that we can reproduce results

##### preprocess data ######
n = dim(banknote.training)[1]  # Size of the training set
m = dim(banknote.test)[1]      # Size of the test set
x = rbind(as.matrix(banknote.training), as.matrix(banknote.test))   
#View(x) ## 200
KK      = 2  ## number of component
p       = dim(x)[2]   ## number of feature
ytrue = c(banknote.training.labels,banknote.test.labels)



##### inspect #####
par(mfrow=c(1,1))
par(mar=c(2,2,2,2)+0.1)
colscale = c("black","red")
pairs(banknote.training,
      col=colscale[banknote.training.labels], 
      pch=as.vector(banknote.training.labels)
      )


####### Initialize the parameters ###########
group1 = banknote.training[banknote.training.labels =="counterfeit",]
group2 = banknote.training[banknote.training.labels =="genuine",]
w          = rep(1,KK)/KK  #Assign equal weight to each component to start with
## eBayes initialize
mu      = array(0, dim=c(KK,p)) 
mu[1,]         = rmvnorm(1, apply(group1,2,mean), cov(group1)) 
mu[2,]         = rmvnorm(1, apply(group2,2,mean), cov(group2)) 
Sigma      = array(0, dim=c(KK,p,p))  #Initial variances are assumed to be the same
Sigma[1,,] = cov(group1)
Sigma[2,,] = cov(group1)
test_cc         = sample(1:KK, m, replace=TRUE, prob=w)
cc = c(as.numeric(banknote.training.labels),test_cc)
cc
# Priors
aa = rep(1, KK)         ## prior of w (follows dirichlet distribution)
dd = apply(x,2,mean)    ## the prior parameters for normal means (in this case,we use the strategy of eBayes)
DD = 10*var(x)          ## the prior variance for normal means 
nu = p                  ## the prior parameters of normal variance (follow inverse wishart distribution(which is a generazation of inverse-gamma))
SS = var(x)/2


# Number of iteration of the sampler
rrr = 5000
burn=1000
# Storing the samples
cc.out    = array(0, dim=c(rrr, n+m))
w.out     = array(0, dim=c(rrr, KK))
mu.out    = array(0, dim=c(rrr, KK, p))
Sigma.out = array(0, dim=c(rrr, KK, p, p))
logpost   = rep(0, rrr)



######## step3: MCMC iterations ##########
for(s in 1:rrr){
  # Sample the indicators in test set
  for(i in (n+1):(n+m)){
    v = rep(0,KK)
    for(k in 1:KK){
      v[k] = log(w[k]) + dmvnorm(x[i,], mu[k,], Sigma[k,,], log=TRUE)  #Compute the log of the weights
    }
    v = exp(v - max(v))/sum(exp(v - max(v)))
    cc[i] = sample(1:KK, 1, replace=TRUE, prob=v)
  }  
  # Sample the weights
  #w=rep(0,KK)
  #w[1] = rbeta(1,aa[1] + sum(cc==1),aa[2] +  sum(cc==2))
  #w[2] = 1-w[1]
  w = as.vector(rdirichlet(1, aa + tabulate(cc)))
  
  # Sample the means
  # reference: lecture P34/73
  # 1) solve: 
  #     function require 2 parameters: solve(a,b)
  #     If missing b, b is taken to be an identity matrix by default,
  #     and solve will return the inverse of a.
  # 2) dd/DD?  
  #    we assume mu follows mu~N(dd,DD)
  #    dd: the prior means for normal means 
  #    DD: the prior variance for normal means 
  DD.st = matrix(0, nrow=p, ncol=p)
  for(k in 1:KK){
    mk    = sum(cc==k)
    xsumk = apply(x[cc==k,], 2, sum)
    DD.st = solve(mk*solve(Sigma[k,,]) + solve(DD))
    dd.st = DD.st%*%(solve(Sigma[k,,])%*%xsumk + solve(DD)%*%dd)
    mu[k,] = as.vector(rmvnorm(1,dd.st,DD.st))
  }
  
  
  # Sample the variances
  # reference: lecture P34/73
  # 1) u/SS
  #    we assume mu follows Sigma~InvWish(u,SS)
  # 2) xcensumk: dim=p*p
  #    xcensumk is the zero centered sum of k-th component.
  #    xcensumk = sum_{i:ci=k} (xi - mu_k)(xi - mu_k)T
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
  for(i in 1:n+m){
    logpost[s] = logpost[s] + log(w[cc[i]]) + mvtnorm::dmvnorm(x[i,],
                                                               mu[cc[i],], 
                                                               Sigma[cc[i],,], 
                                                               log=TRUE)
  }
  
  logpost[s] = logpost[s] + dbeta(w[1],aa[1],aa[2],log=T)
  for(k in 1:KK){
    logpost[s] = logpost[s] + mvtnorm::dmvnorm(mu[k,], dd, DD,log=T)
    logpost[s] = logpost[s] + log(diwish(Sigma[k,,], nu, SS))
  }
  
  if(s/250==floor(s/250)){
    print(paste("s = ", s))
  }
  
}


## Predicted labels of test set
getmode<- function(x){
  return(as.numeric(names(table(x))[table(x) == max(table(x))]))
}

## Plot the logposterior distribution for various samples
par(mfrow=c(1,1))
par(mar=c(4,4,1,1)+0.1)
plot(logpost, type="l", xlab="Iterations", ylab="Log posterior")


## training set

cc.out.df = as.data.frame(cc.out[burn:rrr,])
sum(cc[1:n]!=as.vector(apply(cc.out.df[200:rrr,1:30],2,getmode))) ## 0

## posterior estimate


mu1.out.df = as.data.frame(mu.out[burn:rrr,1,])
mu2.out.df = as.data.frame(mu.out[burn:rrr,2,])
Sigma1.out.df = as.data.frame(Sigma.out[burn:rrr,1,,])
Sigma2.out.df = as.data.frame(Sigma.out[burn:rrr,2,,])

mu1.post = apply(mu1.out.df,2,mean)
mu2.post = apply(mu2.out.df,2,mean)

Sigma1.post = matrix(apply(Sigma1.out.df,2,mean),p,p)
Sigma2.post = matrix(apply(Sigma2.out.df,2,mean),p,p)


mu1.post
mu2.post
## test set
likelihood.group1= dmvnorm(banknote.test,mu1.post,Sigma1.post)
likelihood.group2= dmvnorm(banknote.test,mu2.post,Sigma2.post)
assigned =c()
for (i in 1:length(likelihood.group1)){
  if (likelihood.group1[i] > likelihood.group2[i]){
    assigned = append(assigned,1)
  }else{
    assigned = append(assigned,2)
  }
}
assigned
sum(assigned != ytrue[31:200])
sum(ytrue[30:200]!=as.vector(apply(cc.out.df[500:rrr,30:200],2,getmode)))


# qda
modqda = MASS::qda(grouping=banknote.training.labels, 
                   x=banknote.training, method="mle")
ccpredqda = predict(modqda,newdata=banknote.test)
sum(!(ccpredqda$class == banknote.test.labels)) # 3 error

