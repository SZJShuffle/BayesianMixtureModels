###### Setup data
x = faithful$eruptions
n = length(x)

###### EM algorithm to fit the location-and-scale mixture of 2 Gaussians
w      = 0.5          
mu     = c(mean(x)-sd(x), mean(x)+sd(x))
sigma  = rep(sd(x)/4,2)

s  = 0
sw = FALSE
QQ = -Inf
QQ.out = NULL
epsilon = 10^(-5)

while(!sw){
  ## E step
  v = array(0, dim=c(n,2))
  for(i in 1:n){
    v[i,1] = log(w) + dnorm(x[i], mu[1], sigma[1], log=T)
    v[i,2] = log(1-w) + dnorm(x[i], mu[2], sigma[2], log=T)
    v[i,]  = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,]))) 
  }
  
  ## M step
  # Weights
  w  = mean(v[,1])
  # Means
  mu = rep(0, 2)
  for(k in 1:2){
    for(i in 1:n){
      mu[k]    = mu[k] + v[i,k]*x[i]
    }
    mu[k] = mu[k]/sum(v[,k])
  }
  # Variances
  sigma = rep(0,2)
  for(k in 1:2){
    for(i in 1:n){
      sigma[k] = sigma[k] + v[i,k]*(x[i] - mu[k])^2
    }
    sigma[k] = sqrt(sigma[k]/sum(v[,k]))
  }
  ## Check convergence
  QQn = 0
  for(i in 1:n){
    QQn = QQn + v[i,1]*(log(w) + dnorm(x[i],mu[1],sigma[1],log=TRUE)) +
      v[i,2]*(log(1-w) + dnorm(x[i],mu[2],sigma[2],log=TRUE))
  }
  if(abs(QQn-QQ)/abs(QQn)<epsilon){
    sw=TRUE
  }
  QQ = QQn
  QQ.out = c(QQ.out, QQ)
  s = s + 1
  print(paste(s, QQn))
}


### plot
xx  = seq(0,7,length=150)
nxx = length(xx)
density.EM = rep(0, nxx)
for(s in 1:nxx){
  density.EM[s] = density.EM[s] + w*dnorm(xx[s], mu[1], sigma[1]) + 
    (1-w)*dnorm(xx[s], mu[2], sigma[2])
}
plot(xx, density.EM, col="red", lwd=2, type="l", xlab="Eruptions")
points(x,rep(0,n))


### 
#> w
#[1] 0.348421
#> 1-w
#[1] 0.651579
# > mu
mu ## [1] 2.018646 4.273380
