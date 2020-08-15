### dnorm gives the density, 
### pnorm gives the distribution function, 
### qnorm gives the quantile function, 
### and rnorm generates random deviates.
dnorm(0,0,1) ## PDF,0.398,f(X=0) = 0.398
pnorm(0,0,1) ## CDF,0.5, which means P(X<=0) = 0.5 


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