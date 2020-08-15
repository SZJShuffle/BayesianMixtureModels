
######### 3.Sampling from  hierarchical representation #################
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
