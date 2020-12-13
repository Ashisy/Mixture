## generate data  ##
mydata <-  function(n,lamda,mu,sig)
{
  X1 <-  rnorm(n*lamda[1],mu[1],sig[1]) 
  X2 <-  rnorm(n*lamda[2],mu[2],sig[2])
  dat <-  sample(c(X1,X2)) 
  return(dat)
}
## EM  algorithm  ##
EMnormalmix <-  function(Dat,lamda0,mu0, sig0,
                         max.iter=100, eps=1e-8) 
{ 
  diff=1 
  iter=0
  
  
  while(diff>eps & iter<max.iter) 
  { 
    ## E step ##
    X <- Dat
    
    lamda1 <-  lamda0[1] 
    lamda2 <-  lamda0[2]
    
    mu1 <-  mu0[1]
    mu2 <-  mu0[2]
    
    sig1 <-  sig0[1] 
    sig2 <-  sig0[2]
    
    tmp = lamda1*dnorm(X,mean=mu1,sd=sig1)+ lamda2*dnorm(X,mean=mu2,sd=sig2)
    w1 = lamda1*dnorm(X,mean=mu1,sd=sig1)/tmp 
    w2 = lamda2*dnorm(X,mean=mu2,sd=sig2)/tmp
    
    ## update iteration 
    lamda1.new= mean(w1) 
    lamda2.new= mean(w2) 
    lamda.new = c(lamda1.new, lamda2.new)
    
    ## M step ## 
    mu1.new = sum(w1*X)/sum(w1)
    mu2.new = sum(w2*X)/sum(w2)
    mu.new = c(mu1.new,mu2.new)
    
    sig1square.new = sum(w1*(X-mu1.new)^2)/sum(w1)
    sig2suqare.new = sum(w2*(X-mu2.new)^2)/sum(w2)
    sig.new = c(sqrt(sig1square.new),sqrt(sig2suqare.new))
    
    diff = max(abs(lamda.new-lamda0), abs(mu.new-mu0), abs(sig.new-sig0))
    
    lamda0 = lamda.new 
    mu0 = mu.new 
    sig0 = sig.new
    
    iter=iter+1
  }
  return(list(lamdahat = lamda.new, muhat = mu.new, sighat = sig.new,
              iters = iter,diffs =diff))
}



## simulation ##
set.seed(1) 
n <- 500 
lamda <-  c(0.4,0.6) 
mus <- c(10,15)
sigs <- c(1,3)

## generate data 
D <-  mydata(n,lamda,mus,sigs)

## fit by proposed method 
myfit <- EMnormalmix(D,c(0.1,0.9),c(10,30),c(5,15))
myfit


