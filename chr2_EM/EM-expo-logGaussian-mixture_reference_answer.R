##
w      = 0.1
mu     = mean(log(x))
tau    = sd(log(x))
lambda = 20/mean(x)
## E step
v = array(0, dim=c(n,2))
v[,1] = log(w) + dexp(x, lambda, log=TRUE)    
v[,2] = log(1-w) + dlnorm(x, mu, tau, log=TRUE)
for(i in 1:n){
  v[i,] = exp(v[i,] - max(v[i,]))/sum(exp(v[i,] - max(v[i,])))
}

## M step
w      = mean(v[,1])
lambda = sum(v[,1])/sum(v[,1]*x)
mu     = sum(v[,2]*log(x))/sum(v[,2])
tau    = sqrt(sum(v[,2]*(log(x) - mu)^2)/sum(v[,2]))


## convergence
QQn = 0
for(i in 1:n){
  QQn = QQn + v[i,1]*(log(w) + dexp(x[i], lambda, log=TRUE)) +
    v[i,2]*(log(1-w) + dlnorm(x[i], mu, tau, log=TRUE)) 
}
if(abs(QQn-QQ)/abs(QQn)<epsilon){
  sw=TRUE
}
QQ = QQn

##
w=0.09 lamdba=3.05 mu=0.78 tau=0.38