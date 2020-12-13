## generate data  ##
mydata = function(n,lamda,beta,sig){
  X1 = rnorm(n,2,0.5) 
  X = cbind(1,X1) 
  e1 = rnorm(n,0,sig[1]) 
  e2 = rnorm(n,0,sig[2]) 
  z=rbinom(n,1,lamda[1]) 
  Y=ifelse(z==1,as.numeric(X%*%beta[,1])+e1,as.numeric(X%*%beta[,2])+e2) 
  d = cbind(Y,X)       
  return(d)
}

## EM  algorithm  ##
EMregmix <-  function(data,lamda,beta, sig,
                      max.iter=100, eps=1e-8) 
{ 
  diff=1 
  iter=0
  
  
  while(diff>eps & iter<max.iter) 
  { 
    ## E step ##
    
    Y = data[,1] 
    X = data[,2:3]
    
    
    lamda1 <-  lamda[1] 
    lamda2 <-  lamda[2]
    
    beta1 <- as.matrix(beta[,1]) 
    beta2 <- as.matrix(beta[,2])
    
    sig1 <-  sig[1] 
    sig2 <-  sig[2]
    
    tmp = lamda1*dnorm(Y,mean=as.numeric(X%*%beta1),sd=sig1) +
      lamda2*dnorm(Y,mean=as.numeric(X%*%beta2),sd=sig2)
    
    w1 = lamda1*dnorm(Y,mean=as.numeric(X%*%beta1),sd=sig1)/tmp 
    w2 = lamda2*dnorm(Y,mean=as.numeric(X%*%beta2),sd=sig2)/tmp
    
    W1 = diag(w1) 
    W2 = diag(w2)
    
    ## update iteration 
    lamda1.new=mean(w1) 
    lamda2.new=mean(w2) 
    lamda.new = c(lamda1.new, lamda2.new)
    
    
    ## M step ## 
    beta1.new = as.numeric(solve(t(X)%*%W1%*%X)%*%t(X)%*%W1%*%Y) 
    beta2.new = as.numeric(solve(t(X)%*%W2%*%X)%*%t(X)%*%W2%*%Y) 
    beta.new = c(beta1.new,beta2.new)
    
    
    sig1square.new = t(sqrt(W1)%*%(Y-X%*%beta1.new))%*%
      (sqrt(W1)%*%(Y-X%*%beta1.new))/sum(diag(W1))
    
    sig2suqare.new = t(sqrt(W2)%*%(Y-X%*%beta2.new))%*%
      (sqrt(W2)%*%(Y-X%*%beta2.new))/sum(diag(W2))
    
    sig.new = c(sqrt(sig1square.new),sqrt(sig2suqare.new))
    
    diff = max(abs(lamda.new-lamda), abs(beta.new-beta), abs(sig.new-sig))
    
    lamda = lamda.new 
    beta = beta.new 
    sig = sig.new
    
    
    iter=iter+1
  }
  return(list(lamdahat = lamda.new, betahat = beta.new, sighat = sig.new,
              iters = iter,diffs =diff))
}

## simulation ##
set.seed(1) 
n <- 500 
lamdas <-  c(0.4,0.6) 
betas <- c(10,15)
sigs <- c(1,3)

## generate data 
D <-  mydata(n,lamda,mus,sigs)

## fit by proposed method 
myfit <- EMregmix(D,c(0.1,0.9),c(10,30),c(5,15))
myfit