## stochastic EM algorithm ##
VAREM = function(Dat,lamda0,beta0,tau,p1,max.iter=500) 
{ 
  ## initial ##
  iter=0
  tau = tau
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
    Y = d[,1] 
    X = d[,2:3]
    YY = Dat[,1] 
    XX = Dat[,2:3]
    ## update iteration ##
    
    lamda1.new=mean(p1) 
    lamda2.new=mean(p2)
    
    lamda.new = c(lamda1.new, lamda2.new)
    
    beta1.new <- as.numeric(expTLlaws(X, Y, tau,W0 = p1)$beta)
    beta2.new <- as.numeric(expTLlaws(X, Y, tau,W0 = p2)$beta)
    exVAR1 <- function(d, tau, b)
    {
      d <- data.frame(d)
      d <- subset(d,Z == 1)
      d <- as.matrix(d)
      Y = d[,1] 
      X = d[,2:3]
      ypred <- c(X %*% b)
      w <- as.vector(tau *(Y>= ypred) + (1-tau)* (Y<ypred))
      
      
      
      ## the variance
      Am <- t(X) %*% diag(w) %*% X
      if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}
      
      A <- solve(Am)
      H <- w * diag(X %*% A %*% t(X))
      B <- t(X) %*% diag(w^2*(Y-ypred)^2/(1-H)) %*% X
      Vb <- A %*% B %*% A  # variance function
      return(diag(Vb))
    }
    beta1var <- exVAR1(d,tau,b = beta1.new)
    exVAR2 <- function(d, tau, b)
    {
      d <- data.frame(d)
      d <- subset(d,Z == 0)
      d <- as.matrix(d)
      Y = d[,1] 
      X = d[,2:3]
      ypred <- c(X %*% b)
      w <- as.vector(tau *(Y>= ypred) + (1-tau)* (Y<ypred))
      
      
      
      ## the variance
      Am <- t(X) %*% diag(w) %*% X
      if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}
      
      A <- solve(Am)
      H <- w * diag(X %*% A %*% t(X))
      B <- t(X) %*% diag(w^2*(Y-ypred)^2/(1-H)) %*% X
      Vb <- A %*% B %*% A  # variance function
      return(diag(Vb))
    }
    beta2var <- exVAR2(d,tau,b = beta2.new)
    beta.new = cbind(beta1.new,beta2.new)
    
    ## kernel functiom
    
    h <-  0.05
    kerfun <- function(u){
      1/(h*sqrt(2*pi))*exp(-1/2*(u/h)^2)
    }
    
    
    e1 <- as.numeric(Y-X%*%beta1.new)
    e2 <- as.numeric(Y-X%*%beta2.new)
    i1 <- as.numeric(YY-XX%*%beta1.new)
    i2 <- as.numeric(YY-XX%*%beta2.new)
    I1 <- function(u){u <= 0}
    I2 <- function(u){u > 0}
    t.ini <- 0
    g1.new <- rep(0,length(e1)) 
    
    ## update g1 ##
    v1 <- pnorm(0,e1,h)
    
    a <- sum(p1*I1(e1))
    b <- sum(p1*I2(e1))
    c <- sum(p1*v1*I1(e1))
    dd <- sum(p1*v1*I2(e1))
    bcad <- (b*c-a*dd)
    if (bcad ==0){bcad = 1e-8}
    w11 <-(1-b*(c-a*tau)/bcad)/a
    w21 <- (c-a*tau)/(bcad)
    
    g1.new <- as.numeric(kerfun(outer(i1,e1,"-"))%*%(w11*p1*I1(e1)) + 
                           kerfun(outer(i1,e1,"-"))%*%(w21*p1*I2(e1)))
    
    
    ## update g2 ##
    v2 <- pnorm(0,e2,h)
    
    
    a <- sum(p2*I1(e2))
    b <- sum(p2*I2(e2))
    c <- sum(p2*v2*I1(e2))
    dd <- sum(p2*v2*I2(e2))
    bcad <- (b*c-a*dd)
    if (bcad ==0){bcad = 1e-8}
    w12 <-(1-b*(c-a*tau)/bcad)/a
    w22 <- (c-a*tau)/(bcad)
    
    g2.new <- as.numeric(kerfun(outer(i2,e2,"-"))%*%(w12*p2*I1(e2)) + 
                           kerfun(outer(i2,e2,"-"))%*%(w22*p2*I2(e2)))
    
    
    g.new <-cbind(g1.new,g2.new)
    
    
    lamda0 = lamda.new 
    beta0 = beta.new 
    g0 = g.new
    lamda10 = lamda0[1] 
    lamda20 = lamda0[2]

    iter=iter+1
    
    V_lamda1[iter] <- lamda0[1]
    V_beta11[iter] <- beta1.new[1]
    V_beta12[iter] <- beta1.new[2]
    V_beta21[iter] <- beta2.new[1]
    V_beta22[iter] <- beta2.new[2]
    VV_beta11[iter] <- beta1var[1]
    VV_beta12[iter] <- beta1var[2]
    VV_beta21[iter] <- beta2var[1]
    VV_beta22[iter] <- beta2var[2]
    
    
    g10 = g0[,1] 
    g20 = g0[,2]
    
    
    tmp = lamda10*g10+lamda20*g20
    tmp[tmp == 0] <- 1e-8
    
    p1 = (lamda10*g10)/tmp
    p2 = (lamda20*g20)/tmp
    
    
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

## stochastic EM algorithm ##
VAREM = function(Dat,lamda0,beta0,tau,p0,max.iter=500) 
{ 
  ## initial ##
  iter=0
  tau = tau
  Dat <- Dat
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
    Z <- rMN(p0)
    d <- cbind(Dat,Z)
    
    ## M step ## 
    
    d <-d[sample(nrow(d), nrow(d),replace = T), ]
    p1 <- d[,4]
    p2 <- 1-p1
    Y = d[,1] 
    X = d[,2:3]
    
    ## update iteration ##
    
    lamda1.new=mean(p1) 
    lamda2.new=mean(p2)
    
    lamda.new = c(lamda1.new, lamda2.new)
    
    beta1.new <- as.numeric(expTLlaws(X, Y, tau,W0 = p1)$beta)
    beta2.new <- as.numeric(expTLlaws(X, Y, tau,W0 = p2)$beta)
    exVAR1 <- function(d, tau, b)
    {
      d <- data.frame(d)
      d <- subset(d,Z == 1)
      d <- as.matrix(d)
      Y = d[,1] 
      X = d[,2:3]
      ypred <- c(X %*% b)
      w <- as.vector(tau *(Y>= ypred) + (1-tau)* (Y<ypred))
      
      
      
      ## the variance
      Am <- t(X) %*% diag(w) %*% X
      if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}
      
      A <- solve(Am)
      H <- w * diag(X %*% A %*% t(X))
      B <- t(X) %*% diag(w^2*(Y-ypred)^2/(1-H)) %*% X
      Vb <- A %*% B %*% A  # variance function
      return(diag(Vb))
    }
    beta1var <- exVAR1(d,tau,b = beta1.new)
    exVAR2 <- function(d, tau, b)
    {
      d <- data.frame(d)
      d <- subset(d,Z == 0)
      d <- as.matrix(d)
      Y = d[,1] 
      X = d[,2:3]
      ypred <- c(X %*% b)
      w <- as.vector(tau *(Y>= ypred) + (1-tau)* (Y<ypred))
      
      
      
      ## the variance
      Am <- t(X) %*% diag(w) %*% X
      if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}
      
      A <- solve(Am)
      H <- w * diag(X %*% A %*% t(X))
      B <- t(X) %*% diag(w^2*(Y-ypred)^2/(1-H)) %*% X
      Vb <- A %*% B %*% A  # variance function
      return(diag(Vb))
    }
    beta2var <- exVAR2(d,tau,b = beta2.new)
    beta.new = cbind(beta1.new,beta2.new)
    
    lamda0 = lamda.new 
    beta0 = beta.new 
    iter=iter+1
    
    V_lamda1[iter] <- lamda0[1]
    V_beta11[iter] <- beta1.new[1]
    V_beta12[iter] <- beta1.new[2]
    V_beta21[iter] <- beta2.new[1]
    V_beta22[iter] <- beta2.new[2]
    VV_beta11[iter] <- beta1var[1]
    VV_beta12[iter] <- beta1var[2]
    VV_beta21[iter] <- beta2var[1]
    VV_beta22[iter] <- beta2var[2]
    
    
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

Dat = D
lamda0 = theta0[1:2]
beta0 = matrix(theta0[3:6],2,2)
tau = tau0
p1 = fit$p1
max.iter=500