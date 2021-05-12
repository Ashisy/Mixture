rm(list = ls())

## generate data  ##
mydata = function(n,lamda,beta,sig){
  X1 = runif(n,0,1) 
  X = cbind(1,X1) 
  e1 = rnorm(n,0,sig[1]) 
  e2 = rnorm(n,0,sig[2]) 
  z=c(rep(1,n*lamda[1]),rep(0,n*lamda[2]) )
  Y=ifelse(z==0,as.numeric(X%*%beta[,1])+e1,as.numeric(X%*%beta[,2])+e2) 
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
    beta.new = cbind(beta1.new,beta2.new)
    
    
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
  return(list(lamdahat = lamda.new, betahat = beta.new, p1 = w1))
}

## stochastic EM algorithm ##
VAREM = function(Dat,lamda0,beta0,p1,max.iter=500) 
{ 
  ## initial ##
  iter=0
  ## Store iteration results ##
  V_lamda1 <- rep(0,max.iter)
  V_lamda2 <- rep(0,max.iter)
  V_beta11 <- rep(0,max.iter)
  V_beta12 <- rep(0,max.iter)
  V_beta21 <- rep(0,max.iter)
  V_beta22 <- rep(0,max.iter)
  ## Store beta variance ##  
  VV_beta11 <- rep(0,max.iter)
  VV_beta12 <- rep(0,max.iter)
  VV_beta21 <- rep(0,max.iter)
  VV_beta22 <- rep(0,max.iter)
  
  
  
  while(iter<max.iter) 
  { 
    ## E step ## 
    
    ## initial values 
    
    
    rMN <- function(d)
    {
      i <- 0
      z <- rep(0,length(d))
      for (p in d) 
      {
        i <- i+1
        z[i] <- rbinom(1,1,p)
      }
      return(z)
    }
    Z <- rMN(p1)
    d <- cbind(Dat,Z)
    
    ## M step ## 
    
    d <-d[sample(nrow(d), nrow(d),replace = T), ]
    p1 <- d[,4]
    p2 <- 1-p1
    W1 = diag(p1) 
    W2 = diag(p2)
    Y = d[,1] 
    X = d[,2:3]
    YY = Dat[,1] 
    XX = Dat[,2:3]
    ## update iteration ##
    
    lamda1.new=mean(p1) 
    lamda2.new=mean(p2)
    
    lamda.new = c(lamda1.new, lamda2.new)
    
    ## M step ## 
    beta1.new = as.numeric(solve(t(X)%*%W1%*%X)%*%t(X)%*%W1%*%Y) 
    beta2.new = as.numeric(solve(t(X)%*%W2%*%X)%*%t(X)%*%W2%*%Y)
    lmVAR1 <- function(dat)
    {
      d <- data.frame(dat)
      d <- subset(d,Z == 1)
      d <- as.matrix(d)
      y <- dat[, 1]
      x <- dat[,3]
      d <- data.frame(cbind(y,x))
      fit <- lm(y~x,d)
      var <- as.numeric(diag(vcov(fit))) 
      return(var)
    }
    beta1var <- lmVAR1(dat = d)
    lmVAR2 <- function(dat)
    {
      d <- data.frame(dat)
      d <- subset(d,Z == 0)
      d <- as.matrix(d)
      y <- dat[, 1]
      x <- dat[,3]
      d <- data.frame(cbind(y,x))
      fit <- lm(y~x,d)
      var <- as.numeric(diag(vcov(fit))) 
      return(var)
    }
    beta2var <- lmVAR2(dat = d)
    
    beta.new = cbind(beta1.new,beta2.new)
    
    
    sig1square.new = t(sqrt(W1)%*%(Y-X%*%beta1.new))%*%
      (sqrt(W1)%*%(Y-X%*%beta1.new))/sum(diag(W1))
    
    sig2suqare.new = t(sqrt(W2)%*%(Y-X%*%beta2.new))%*%
      (sqrt(W2)%*%(Y-X%*%beta2.new))/sum(diag(W2))
    
    sig.new = c(sqrt(sig1square.new),sqrt(sig2suqare.new))
    
    lamda = lamda.new 
    beta = beta.new 
    sig = sig.new
    
    
    iter=iter+1
    
    lamda0 = lamda.new 
    lamda1 = lamda0[1] 
    lamda2 = lamda0[2]
    
    
    V_lamda1[iter] <- lamda0[1]
    V_beta11[iter] <- beta1.new[1]
    V_beta12[iter] <- beta1.new[2]
    V_beta21[iter] <- beta2.new[1]
    V_beta22[iter] <- beta2.new[2]
    VV_beta11[iter] <- beta1var[1]
    VV_beta12[iter] <- beta1var[2]
    VV_beta21[iter] <- beta2var[1]
    VV_beta22[iter] <- beta2var[2]
    
    
   
    beta1 <- as.matrix(beta[,1]) 
    beta2 <- as.matrix(beta[,2])
    
    sig1 <-  sig[1] 
    sig2 <-  sig[2]
    
    tmp = lamda1*dnorm(YY,mean=as.numeric(XX%*%beta1),sd=sig1) +
      lamda2*dnorm(YY,mean=as.numeric(XX%*%beta2),sd=sig2)
    tmp[tmp == 0] <- 1e-8
    p1 = lamda1*dnorm(YY,mean=as.numeric(XX%*%beta1),sd=sig1)/tmp 
    p2 = lamda2*dnorm(YY,mean=as.numeric(XX%*%beta2),sd=sig2)/tmp
  } 
  V_lamda2 <- 1 - V_lamda1
  VAR_lamda1 <- drop(t(V_lamda1-mean(V_lamda1))%*%((V_lamda1-mean(V_lamda1))))/(max.iter-1)*(1 + 1/max.iter)+ 
    drop(t(V_lamda1)%*%V_lamda1)/(max.iter*nrow(Dat))
  
  VAR_lamda2 <- VAR_lamda1
  
  VAR_beta11 <- drop(t(V_beta11-mean(V_beta11))%*%((V_beta11-mean(V_beta11))))/(max.iter-1)*(1 + 1/max.iter) + sum(VV_beta11)/max.iter
  VAR_beta12 <- drop(t(V_beta12-mean(V_beta12))%*%((V_beta12-mean(V_beta12))))/(max.iter-1)*(1 + 1/max.iter) + sum(VV_beta12)/max.iter
  VAR_beta21 <- drop(t(V_beta21-mean(V_beta21))%*%((V_beta21-mean(V_beta21))))/(max.iter-1)*(1 + 1/max.iter) + sum(VV_beta21)/max.iter
  VAR_beta22 <- drop(t(V_beta22-mean(V_beta22))%*%((V_beta22-mean(V_beta22))))/(max.iter-1)*(1 + 1/max.iter) + sum(VV_beta22)/max.iter
  
  return(c(VAR_lamda1,VAR_lamda2,VAR_beta11,VAR_beta12,VAR_beta21,VAR_beta22))
}

## simulation ## 
n <- 500
theta0 <- c(0.5,0.5,10,-10,-10,10) # true value
tau0= 0.5# true value
D <- mydata(n,theta0[1:2],matrix(theta0[3:6],2,2) ,sig = c(1,1))
plot(D[, 3],D[, 1])
theta.ini <- c(0.5,0.5,3,6,8,-10)
fit <- EMregmix(data = D,lamda =theta.ini[1:2],beta =matrix(theta.ini[3:6],2,2),
                sig = c(1,1))
sqrt(VAREM(Dat = D,lamda0 = theta0[1:2],
           beta0 = matrix(theta0[3:6],2,2),p1 = fit$p1))









## simulation
set.seed(123)
n <- 500
theta0 <- c(0.5,0.5,10,-10,-10,10) # true value
tau0 <- 0.5
NS <- 50

p <- length(theta0)
theta.est <- matrix(0, NS, p)
theta.ese <- matrix(0, NS, p)

for (kk in 1:NS){
  ## generated data
  
  D <- mydata(n,theta0[1:2],matrix(theta0[3:6],2,2) ,sig = c(1,1))
  if (kk==1){
    y <- D[, 1]
    x <- D[, 3]
    plot(x,y)
  }
  
  theta.ini <- c(0.5,0.5,3,6,8,10)
  
  fit <- EMregmix(data = D,lamda =theta.ini[1:2],beta =matrix(theta.ini[3:6],2,2),
                  sig = c(1,1))
  variance <- VAREM(Dat = D,lamda0 = theta0[1:2],
                    beta0 = matrix(theta0[3:6],2,2),p1 = fit$p1)
  
  theta.est[kk,1:2 ] <- fit$lamdahat
  theta.est[kk,3:6 ] <- as.vector(fit$betahat)
  theta.ese[kk, ] <- sqrt(variance)
}

### 
bias <- apply(theta.est, 2, mean) - theta0
sd <- apply(theta.est, 2, sd)
ese <- apply(theta.ese, 2, mean)
alp <- 0.05
c.alp <- qnorm(1-alp/2, 0, 1)
theta000 <- matrix( rep(theta0, NS), nrow=NS, ncol=6, byrow = TRUE)
coverage <-  (theta.est-c.alp* theta.ese <= theta000) *
  (theta000 <= theta.est + c.alp* theta.ese )
CP <- apply(coverage, 2, mean)  


out <- rbind(bias, sd, ese, CP)
colnames(out) <- c("lamba1", "lamba2", "beta11","beta12","beta21","beta22")

out

