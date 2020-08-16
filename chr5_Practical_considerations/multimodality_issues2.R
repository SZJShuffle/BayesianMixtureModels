## Illustrating that the EM might fail for numerical reasons if a component is “numerically empty”

### Loading data and setting up global variables
rm(list=ls())
library(mclust)
library(mvtnorm)
library(ellipse)

## Setup data
data(iris)
x       = as.matrix(iris[,-5])
n       = dim(x)[1]
p       = dim(x)[2]       # Number of features
KK      = 3
epsilon = 0.00000001

# Initialize the parameters of the algorithm
w   = rep(1,KK)/KK  #Assign equal weight to each component to start with
mu  = matrix(0, KK, p)  # Initialize in the true values
mu[1,] = apply(x, 2, mean)
mu[2,] = apply(x, 2, mean) + c(2.2, 2.2, 2.2, 2.2)
mu[3,] = apply(x, 2, mean) + c(-2.2, -2.2, -2.2, -2.2)
Sigma      = array(0, dim=c(KK,p,p))  #Initial variances are assumed to be the same
Sigma[1,,] = var(x)/3
Sigma[2,,] = var(x)/3
Sigma[3,,] = var(x)/3

# Plot the data along with the estimates of the components
colscale = c("black","blue","red")
par(mfrow=c(p,p))
for(k in 1:p){
  for(l in 1:p){
    if(k!=l){
      par(mar=c(3,3,1,1)+0.1)
      plot(x[,k], x[,l], type="n", xlab="", ylab="", xlim=c(min(c(x[,k], mu[,k])),max(c(x[,k], mu[,k]))), ylim=c(min(c(x[,l], mu[,l])),max(c(x[,l], mu[,l]))))
      for(r in 1:KK){
        lines(ellipse(x=Sigma[r,c(k,l),c(k,l)], centre=mu[r,c(k,l)], level=0.50), col="gold1", lty=1, lwd=1)
        lines(ellipse(x=Sigma[r,c(k,l),c(k,l)], centre=mu[r,c(k,l)], level=0.82), col="gold1", lty=1, lwd=1)
        lines(ellipse(x=Sigma[r,c(k,l),c(k,l)], centre=mu[r,c(k,l)], level=0.95), col="gold1", lty=1, lwd=1)
      }
      text(x[,k], x[,l], labels=as.numeric(iris[,5]), col=colscale[iris[,5]])
      points(mu[,k], mu[,l], pch=19, col="gold1", cex=2)
    }else{
      plot(seq(0,5), seq(0,5), type="n", xlab="", ylab="", axes=FALSE)
      text(2.5,2.5,colnames(x)[k], cex=1.5)
    }
  }
}

## Run the EM algorithm.  It will fail in the first iteration
sw     = FALSE
QQ     = -Inf
QQ.out = NULL
s      = 0

while(!sw){
  ## E step
  v = array(0, dim=c(n,KK))
  for(k in 1:KK){  #Compute the log of the weights
    v[,k] = log(w[k]) + dmvnorm(x, mu[k,], Sigma[k,,], log=TRUE) 
  }
  for(i in 1:n){
    v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))  #Go from logs to actual weights in a numerically stable manner
  }
  
  ## M step
  w = apply(v,2,mean)
  mu = matrix(0, nrow=KK, ncol=p)
  for(k in 1:KK){
    for(i in 1:n){
      mu[k,]    = mu[k,] + v[i,k]*x[i,]
    }
    mu[k,] = mu[k,]/sum(v[,k])
  }
  Sigma = array(0,dim=c(KK, p, p))
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
}

QQn
