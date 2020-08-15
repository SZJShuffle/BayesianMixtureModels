## Example of an EM algorithm for fitting a mixture of Exponential and log-Gaussian mixture model 

## Clear the environment and load required libraries
rm(list=ls())
set.seed(81196)    # So that results are reproducible (same simulated data every time)

## 1.load data 
fuses = read.csv('fuses.csv',stringsAsFactors = F,header = F)
X = fuses$V1

## 3.Plot the Empirical distribution
par(mfrow=c(1,1))
hist(X,breaks=50)


## 4.Initialize the parameters
KK=2            #number of cluster
w = 0.05       
lambda    = 0.1   #param of Exponential-distribution.
mu = mean(log(X))    
tau = sd(X) 
n=length(X)


## Plot the initial guess for the density
xx = seq(0,11,length=200)
yy = w*dexp(xx,rate = lambda) + (1-w)*dlnorm(xx,meanlog = mu,sdlog = sigma)
plot(xx, yy, type="l", ylim=c(0, max(yy)), xlab="x", ylab="Initial density")


s  = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)


##Checking convergence of the algorithm
while(!sw){
  ## E step
  v = array(0, dim=c(n,KK))
  
  ## log form
  #v[,1] = log(w)   + dexp(X, rate=lambda, log=TRUE)     #Compute the log of the weights
  #v[,2] = log(1-w) + dlnorm(X, mean=mu, sd=tau, log=TRUE)  #Compute the log of the weights
  #for(i in 1:n){
  #  v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))  #Go from logs to actual weights in a numerically stable manner
  #}
  
  ## raw form (both OK)
  v[, 1] = w * dexp(X, lambda)
  v[, 2] = (1-w) * dlnorm(X, mu, tau)
  v = v / (v[, 1] + v[, 2])
  
  ## M step
  # Weights
  w = mean(v[,1])
  lambda=sum(v[,1])/sum(v[,1]*X)
  mu = sum(v[,2]*log(X))/sum(v[,2])
  tau = sqrt(sum(v[,2]*(log(X)-mu)^2)/sum(v[,2]))
  
  ##Check convergence
  QQn = 0
  for(i in 1:n){
    QQn = QQn + v[i,1]*(log(w) + dexp(X[i],rate = lambda,log=T)) +
      v[i,2]*(log(1-w) + dlnorm(X[i], mu, tau, log=T))
  }
  if(abs(QQn-QQ)/abs(QQn)<epsilon){
    sw=TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))
  
  #Plot current estimate over data
  layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
  par(mar=c(3.1,4.1,0.5,0.5))
  plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)
  
  ## 
  par(mar=c(5,4,1.5,0.5))
  xx = seq(0,11,length=200)
  yy = w*dexp(xx,rate=lambda) + (1-w)*dlnorm(xx,meanlog = mu,sdlog = tau)
  plot(xx,yy,type='l',xlab="x", ylab="PDF")  
}


#Plot current estimate over data
layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
par(mar=c(3.1,4.1,0.5,0.5))
plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)

## 
par(mar=c(5,4,1.5,0.5))
xx = seq(0,11,length=200)
yy = w*dexp(xx,rate=lambda) + (1-w)*dlnorm(xx,meanlog = mu,sdlog = sigma)
plot(xx,yy,type='l',xlab="x", ylab="PDF")  

round(w,2)
round(mu,2)
round(tau,2)
round(lambda,2)
