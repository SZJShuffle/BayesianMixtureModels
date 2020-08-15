
## guess
## E step
w      = 0.5
lambda = mean(x)

## E step
v = array(0, dim=c(n,2))
for(i in 1:n){
  if(x[i]==0){
    v[i,1] = w    
    v[i,2] = (1-w)*dpois(x[i], lambda)  
    v[i,]  = v[i,]/sum(v[i,])
  }else{
    v[i,1] = 0
    v[i,2] = 1
  }
}


w = mean(v[,1])
lambda = sum(x)/sum(v[,2])

##Check convergence
QQn = 0
for(i in 1:n){
  if(x[i]==0){
    QQn = QQn + v[i,1]*log(w) + v[i,2]*(log(1-w) + dpois(x[i], lambda, log=TRUE))
  }else{
    QQn = QQn + v[i,2]*(log(1-w) + dpois(x[i], lambda, log=TRUE))
  }
}
if(abs(QQn-QQ)/abs(QQn)<epsilon){
  sw=TRUE
}
QQ = QQn

