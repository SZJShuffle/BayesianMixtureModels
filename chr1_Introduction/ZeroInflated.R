
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