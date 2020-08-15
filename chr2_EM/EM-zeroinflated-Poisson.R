#### Example of an EM algorithm for fitting a zero-inflated Poisson 

## Clear the environment and load required libraries
rm(list=ls())
set.seed(81196)    # So that results are reproducible (same simulated data every time)


## 1.load data 
nestsize = read.csv('nestsize.csv',stringsAsFactors = F)
X = nestsize$X4


## 3.Plot the Empirical distribution
table(X)
par(mfrow=c(1,1))
barplot(table(X))

##### Run the actual EM algorithm ######

## 4.Initialize the parameters
KK=2            ## number of cluster
w     = 1/2     #Assign equal weight to each component to start with
lambda    = 3   #Random cluster centers randomly spread over the support of the data
n=length(X)

## 5.Plot the initial guess for the density
xx = seq(0,11)
yy = w*c(1,rep(0,length(xx)-1)) + (1-w)*dpois(xx,lambda = lambda)
barplot(yy, ylim=c(0, max(yy)), xlab="x", ylab="Initial Mass")

s  = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)

## 6.Checking convergence of the algorithm

while(!sw){
  ## E step
  v = array(0, dim=c(n,KK))
  
  # p(zi=1|xi,THETA) => v[,1]
  # p(zi=2|xi,THETA) => v[,1]
  v[,1] = ifelse(X>0,0,w/(w + (1-w)*dpois(X,lambda=lambda)))
  v[,2] = ifelse(X>0,1,(1-w)*exp(-lambda)/(w + (1-w)*dpois(X,lambda=lambda)))
  for(i in 1:n){
    v[i,] = v[i,]/sum(v[i,])
  }
  
  ## M step
  # Weights and lambda
  w = mean(v[,1])
  lambda = sum(X)/sum(v[,2])

  ##Check convergence
  ## QQn:Q fucntion value at current iteration
  ## Q:Q fucntion value at previous iteration
  QQn = 0
  for(i in 1:n){
    QQn = QQn + v[i,1]*(log(w)) +
      v[i,2]*(log(1-w) + dpois(X[i], lambda = lambda, log=TRUE))
  }
  print(QQn)
  print(QQ)
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
 
  par(mar=c(5,4,1.5,0.5))
  xx = seq(0,11)
  yy = w*c(1,rep(0,length(xx)-1)) + (1-w)*dpois(xx,lambda = lambda)
  barplot(yy, ylim=c(0, max(yy)), xlab="x", ylab="PMF")
}

print(paste(round(lambda,digits=2),round(w,digits=2)))

#Plot final estimate over data
layout(matrix(c(1,2),2,1), widths=c(1,1), heights=c(1.3,3))
par(mar=c(3.1,4.1,0.5,0.5))
plot(QQ.out[1:s],type="l", xlim=c(1,max(10,s)), las=1, ylab="Q", lwd=2)

par(mar=c(5,4,1.5,0.5))
xx = seq(0,11)
yy = w*c(1,rep(0,length(xx)-1)) + (1-w)*dpois(xx,lambda = lambda)
barplot(yy, ylim=c(0, max(yy)), xlab="x", ylab="Final PMF")

