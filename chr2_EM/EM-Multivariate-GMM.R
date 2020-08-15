#### Example of an EM algorithm for fitting a mixtures of K p-variate Gaussian components
#### The algorithm is tested using simulated data

## Clear the environment and load required libraries
rm(list=ls())
library(mvtnorm)    # Multivariate normals are not default in R
library(ellipse)    # Required for plotting
set.seed(63252)     # For reproducibility

## Generate data from a mixture with 3 components
KK      = 3
p       = 2
w.true = c(0.5,0.3,0.2)  # True weights associated with the components
mu.true     = array(0, dim=c(KK,p))
mu.true[1,] = c(0,0)   #True mean for the first component
mu.true[2,] = c(5,5)   #True mean for the second component
mu.true[3,] = c(-3,7)   #True mean for the third component
Sigma.true      = array(0, dim=c(KK,p,p))
Sigma.true[1,,] = matrix(c(1,0,0,1),p,p)   #True variance for the first component
Sigma.true[2,,] = matrix(c(2,0.9,0.9,1),p,p)   #True variance for the second component
Sigma.true[3,,] = matrix(c(1,-0.9,-0.9,4),p,p)   #True variance for the third component
n  = 120
cc = sample(1:3, n, replace=T, prob=w.true)
x  = array(0, dim=c(n,p))
for(i in 1:n){
  x[i,] = rmvnorm(1, mu.true[cc[i],], Sigma.true[cc[i],,])
}

par(mfrow=c(1,1))
plot(x[,1], x[,2], col=cc, type="n", xlab=expression(x[1]), ylab=expression(x[2]))
text(x[,1], x[,2], seq(1,n), col=cc, cex=0.6)
for(k in 1:KK){
  lines(ellipse(x=Sigma.true[k,,], centre=mu.true[k,], level=0.50), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma.true[k,,], centre=mu.true[k,], level=0.82), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma.true[k,,], centre=mu.true[k,], level=0.95), col="grey", lty=2, lwd=2)
}
#title(main="Data + True Components")


### Run the EM algorithm
## Initialize the parameters
w   = rep(1,KK)/KK  #Assign equal weight to each component to start with
mu  = rmvnorm(KK, apply(x,2,mean), var(x))   #RandomCluster centers randomly spread over the support of the data
Sigma      = array(0, dim=c(KK,p,p))  #Initial variances are assumed to be the same
Sigma[1,,] = var(x)/KK  
Sigma[2,,] = var(x)/KK
Sigma[3,,] = var(x)/KK

par(mfrow=c(1,1))
plot(x[,1], x[,2], col=cc, xlab=expression(x[1]), ylab=expression(x[2]))
for(k in 1:KK){
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.50), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.82), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.95), col="grey", lty=2, lwd=2)
}
title(main="Initial estimate + Observations")

s       = 0
sw      = FALSE
QQ      = -Inf
QQ.out  = NULL
epsilon = 10^(-5)

while(!sw){
  ## E step
  v = array(0, dim=c(n,KK))
  for(k in 1:KK){
    v[,k] = log(w[k]) + dmvnorm(x, mu[k,], Sigma[k,,],log=TRUE)  #Compute the log of the weights
  }
  for(i in 1:n){
    v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))  #Go from logs to actual weights in a numerically stable manner
  }
  
  ## M step
  w     = apply(v,2,mean)
  mu    = array(0, dim=c(KK, p))
  for(k in 1:KK){
    for(i in 1:n){
      mu[k,]    = mu[k,] + v[i,k]*x[i,]
    }
    mu[k,] = mu[k,]/sum(v[,k])
  }
  Sigma = array(0, dim=c(KK, p, p))
  for(k in 1:KK){
    for(i in 1:n){
      Sigma[k,,] = Sigma[k,,] + v[i,k]*(x[i,] - mu[k,])%*%t(x[i,] - mu[k,])
    }
    Sigma[k,,] = Sigma[k,,]/sum(v[,k])
  }
  
  ##Check convergence
  QQn = 0
  for(i in 1:n){
    for(k in 1:KK){
      QQn = QQn + v[i,k]*(log(w[k]) + dmvnorm(x[i,],mu[k,],Sigma[k,,],log=TRUE))
    }
  }
  if(abs(QQn-QQ)/abs(QQn)<epsilon){
    sw=TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))
  
  #Plot current components over data
  layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
  par(mar=c(3.1,4.1,0.5,0.5))
  plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q")
  
  par(mar=c(5,4,1,0.5))
  plot(x[,1], x[,2], col=cc, main=paste("s =",s,"   Q =", round(QQ.out[s],4)), 
       xlab=expression(x[1]), ylab=expression(x[2]), lwd=2)
  for(k in 1:KK){
    lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.50), col="grey", lty=2, lwd=2)
    lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.82), col="grey", lty=2, lwd=2)
    lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.95), col="grey", lty=2, lwd=2)
  }
}

#Plot current components over data
layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
par(mar=c(3.1,4.1,0.5,0.5))
plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)

par(mar=c(5,4,1,0.5))
plot(x[,1], x[,2], col=cc, main=paste("s =",s,"   Q =", round(QQ.out[s],4)), xlab=expression(x[1]), ylab=expression(x[2]))
for(k in 1:KK){
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.50), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.82), col="grey", lty=2, lwd=2)
  lines(ellipse(x=Sigma[k,,], centre=mu[k,], level=0.95), col="grey", lty=2, lwd=2)
}
