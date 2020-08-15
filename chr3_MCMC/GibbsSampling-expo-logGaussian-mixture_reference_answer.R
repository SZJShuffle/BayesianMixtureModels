w      = 0.1
mu     = mean(log(x))
tau    = sd(log(x))
lambda = 20/mean(x)
cc     = sample(1:2, n, TRUE, c(1/2, 1/2))


# Full conditional for cc
v = rep(0,2)
for(i in 1:n){
  v[1]  = log(w) + dexp(x[i], lambda, log=TRUE)
  v[2]  = log(1-w) + dlnorm(x[i], mu, tau, log=TRUE)
  v     = exp(v - max(v))/sum(exp(v - max(v)))
  cc[i] = sample(1:2, 1, replace=TRUE, prob=v)
}

# Full conditional for w
w = rbeta(1, 1+sum(cc==1), 1+n-sum(cc==1))
# Full conditional for w
w = rbeta(1, 1+sum(cc==1), 1+sum(cc==2))

# Full conditional for lambda
lambda = rgamma(1, 1 + sum(cc==1), 1 + sum(x[cc==1]))

# Full conditional for mu
mean.post = (sum(log(x[cc==2]))/tau^2 + 0)/(sum(cc==2)/tau^2 + 1)
std.post = sqrt(1/(sum(cc==2)/tau^2 + 1))
mu = rnorm(1, mean.post, std.post)

# Full conditional for tau
tau = sqrt(1/rgamma(1, 2 + sum(cc==2), 1 + sum((log(x[cc==2]) - mu)^2)))
                
                