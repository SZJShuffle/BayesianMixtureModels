n        = length(x)
cc       = rep(0, n)
cc[x==0] = sample(1:2, sum(x==0), replace=T, prob=c(1/2, 1/2))
cc[x!=0] = 2
lambda   = mean(x)
w        = 0.2


# Full conditional for cc
for(i in 1:n){
  v = rep(0,2)
  if(x[i]==0){
    v[1] = log(w)
    v[2] = log(1-w) + dpois(x[i], lambda, log=TRUE)
    v    = exp(v - max(v))/sum(exp(v - max(v)))
  }else{
    v[1] = 0
    v[2] = 1
  }
  cc[i] = sample(1:2, 1, replace=TRUE, prob=v)
}

# Full conditional for w
w = rbeta(1, 1+sum(cc==1), 1+sum(cc==2))

lambda = rgamma(1, sum(x[cc==2]) + 1, sum(cc==2) + 1)
