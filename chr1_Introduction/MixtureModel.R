
### dnorm gives the density, 
### pnorm gives the distribution function, 
### qnorm gives the quantile function, 
### and rnorm generates random deviates.
dnorm(0,0,1) ## PDF,0.398,f(X=0) = 0.398
pnorm(0,0,1) ##CDF,0.5,表示P(X<=0) = 0.5 


############ 1.GMM #################
## mixture of 2 normal dist
## 1.1multi-modal:large difference in mean, little difference in var
x = seq(-5, 12, length=100)
y = 0.6*dnorm(x, 0, 1) + 0.4*dnorm(x, 5, 2)
par(mar=c(4,4,1,1)+0.1)
plot(x, y, type="l", ylab="Density", las=1, lwd=2)


## 1.2skewed: Mixture of univariate Gaussians, unimodal skewed
x = seq(-5, 12, length=100)
y = 0.55*dnorm(x, 0, sqrt(2)) + 0.45*dnorm(x, 3, 4)
par(mar=c(4,4,1,1)+0.1)
plot(x, y, type="l", ylab="Density", las=1, lwd=2)


## 1.3unimodal: heavy tail
x = seq(-12, 12, length=100)
## y和z的期望与方差都一样，但是y是单峰肥尾分布
y = 0.4*dnorm(x, 0, sqrt(2)) + 0.4*dnorm(x, 0,sqrt(16)) +  0.2*dnorm(x, 0, sqrt(20))
z = dnorm(x, 0, 0.4*sqrt(2) + 0.4*sqrt(16) + 0.2*sqrt(20))
par(mar=c(4,4,1,1)+0.1)
plot(x, y, type="l", ylab="Density", las=1, lwd=2)
lines(x,z,lty=2)


############ 2.zero-inflated #################
## Zero inflated negative binomial distribution
x = seq(0, 15)
y = dnbinom(x, 8, 0.6)
z = 0.2*c(1,rep(0,length(x)-1)) + (1-0.2)*y
par(mfrow=c(2,1))
par(mar=c(4,4,2,2)+0.1)
barplot(y, names.arg=x, las=1, xlab = "x", ylab="Probability", 
        border=NA, main="Negative Binomial")
par(mar=c(4,4,1,1)+0.1)
barplot(z, names.arg=x, las=1, xlab = "x", ylab="Probability", 
        border=NA, main="Zero-inflated Negative Binomial")

## Zero-inflated log Gaussian distribution
x = seq(-2, 15, length=1000)
y = plnorm(x, 1.5, 0.5) ## cdf
z = 0.3*as.numeric(x>=0) + (1-0.3)*y
par(mar=c(4,4,1,1)+0.1)
plot(x, y, type="l", las=1, lty=2, xlab="x", 
     ylab="Cumulative distribution Function", lwd=2)
lines(x, z, lty=1, lwd=2)
legend(4, 0.45, c("Zero infla. log Gaussian","log Gaussian"), 
       lty=c(1,2), bty="n", lwd=c(2,2))



## zero-inflated gaussain
x = seq(-15, 15)
y = dnorm(x, 8, 5)
z = 0.2*as.numeric(x ==0) + (1-0.2)*y
par(mfrow=c(2,1))
par(mar=c(4,4,2,2)+0.1)
barplot(y, names.arg=x, las=1, xlab = "x", ylab="Probability", 
        border=NA, main="Gaussian")
par(mar=c(4,4,1,1)+0.1)
barplot(z, names.arg=x, las=1, xlab = "x", ylab="Probability", 
        border=NA, main="Zero-inflated Gaussian")


######### 3，sampling from  hierarchical representation #################
# Generate n observations from a mixture of two Gaussian 
# distributions
n     = 50           # Size of the sample to be generated
w     = c(0.6, 0.4)  # Weights
mu    = c(0, 5)      # Means
sigma = c(1, 2)      # Standard deviations
cc    = sample(1:2, n, replace=T, prob=w)
x     = rnorm(n, mu[cc], sigma[cc])

# Plot f(x) along with the observations 
# just sampled
xx = seq(-5, 12, length=200)
yy = w[1]*dnorm(xx, mu[1], sigma[1]) + 
  w[2]*dnorm(xx, mu[2], sigma[2])
par(mar=c(4,4,1,1)+0.1)
plot(xx, yy, type="l", ylab="Density", xlab="x", las=1, lwd=2)
points(x, y=rep(0,n), pch=1)

# Generate n observations from a mixture of 3 poisson
# sample 200 random numbers from a mixture of 3 Poisson distributions 
# with means 1, 2 and 6 and weights 0.7, 0.2 and 0.1, respectively, 
# and generate a barplot with the empirical frequencies of all the integers 
# included in the sample.
n     = 200           # Size of the sample to be generated
w     = c(0.7, 0.2,0.1)  # Weights
lambda    = c(1, 2,6)      # Means
cc    = sample(1:3, n, replace=T, prob=w)
x     = rpois(n, lambda[cc])

# Plot f(x) along with the observations 
par(mar=c(4,4,1,1)+0.1)
barplot(table(x)/n,ylab="Freq")

