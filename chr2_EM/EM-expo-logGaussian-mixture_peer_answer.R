## Clear the environment and load required libraries
rm(list=ls())

# read data
x = read.csv("fuses.csv", header = FALSE)
x = as.numeric(x[,1])

# explore data
hist(x, breaks=length(x)/10, freq=FALSE)

## Run the actual EM algorithm
## Initialize the parameters
# by visual inspection, it looks like that for x < 1, data is indeed
# very well exponentially distributed, and for x > 1, lognormal distributed
# (or similar distributions).
# Most of the mass in on the lognormal.
w     = length(x[x < 1]) / length(x)
# since lambda is the mean, start with that for lambda
lambda = mean(x[x < 1])
# since the lognormal median is exp(mu), initialize mu accordingly
mu = log(median(x[x >= 1]))
# singe the lognormal mean is exp(mu + tau^2/2), solve for it
tau = ( 2  * ( log(mean(x[x>=1])) - mu ) ) ^(0.5)

# plot current density and data
plotd <- function(s) {
  xx = seq(0,max(x),length=200)
  yy = w * dexp(xx, lambda) + (1-w) * dlnorm(xx, mean=mu, sd=tau)
  hist(x, breaks=length(x)/10, freq=FALSE, main=paste0("s=",s))
  lines(xx, yy, type="l", col="red", lwd=2)
}
plotd(0)

s  = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)

##Checking convergence of the algorithm
n = length(x)
KK = 2
while(!sw){
  ## E step
  v = array(0, dim=c(n,KK))
  v[,1] = log(w)   + dexp(x, rate=lambda, log=TRUE)     #Compute the log of the weights
  v[,2] = log(1-w) + dlnorm(x, mean=mu, sd=tau, log=TRUE)  #Compute the log of the weights
  for(i in 1:n){
    v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))  #Go from logs to actual weights in a numerically stable manner
  }
  
  ## M step
  # Weights
  w = mean(v[,1]) # always this, whatever the pdf!
  
  # note: expressions could be vectorized, keeping the loops
  # to be consistent with original text
  
  # lambda
  numerator = 0
  denominator = 0
  for(i in 1:n){
    numerator   = numerator   + v[i,1]
    denominator = denominator + v[i,1]*x[i]
  }
  lambda = numerator / denominator
  
  # mu
  numerator = 0
  denominator = 0
  for(i in 1:n){
    numerator   = numerator   + v[i,2]*log(x[i])
    denominator = denominator + v[i,2]
  }
  mu = numerator / denominator
  
  # tau
  numerator = 0
  denominator = 0
  for(i in 1:n){
    numerator   = numerator   + v[i,2]*((log(x[i]) - mu)^2)
    denominator = denominator + v[i,2]
  }
  tau = sqrt(numerator/denominator)
  
  ##Check convergence
  QQn = 0
  for(i in 1:n){ # done, ok?
    QQn = QQn + v[i,1]*(log(w)   + dexp(x[i], rate=lambda, log=TRUE)) +
      v[i,2]*(log(1-w) + dlnorm(x[i], mean=mu, sd=tau, log=TRUE))
  }
  if(abs(QQn-QQ)/abs(QQn)<epsilon){
    sw=TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))
  
  #Plot current estimate over data (final for last loop)
  plotd(s)
}

##  "w=0.09 lamdba=3.05 mu=0.78 tau=0.38"
print(paste0("w=", round(w,2), " lamdba=",round(lambda,2), " mu=", round(mu,2), " tau=", round(tau,2)))