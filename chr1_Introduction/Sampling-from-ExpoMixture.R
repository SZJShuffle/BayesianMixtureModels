# Generate n observations from a mixture of 4 exponential distribution 
n     = 100           # Size of the sample to be generated
w     = c(0.3, 0.25,0.25,0.2)  # Weights
lambda    = c(1, 4,7,10)      # Means

cc    = sample(1:4, n, replace=T, prob=w)
x     = rexp(n, 1/lambda[cc])
mean(x) ## 6.07415
var(x) ## 49.15045
