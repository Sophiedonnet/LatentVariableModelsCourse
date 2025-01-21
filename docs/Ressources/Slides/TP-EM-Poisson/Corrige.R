n = 1000
omega = c(1/3,1/4,1-1/3-1/4)
mu = c(5,10,20)



### 1 .  Simulate data

 
Z = sample(1:3,n,replace= TRUE,prob   = omega)
muZ = mu[Z]
y = rpois(n,mu[Z])

### 2. Hitogramm

hist(y,nclass=40)


### 3. Write the function to compute de log likelihood


LogLikPoissonMixture <- function(y,mu,omega){
  K <- length(mu)
  lik <-  rowSums(sapply(1:K, function(k){dpois(y,mu[k])*omega[k]}))
  return(sum(log(lik)))
}


### 4.  Write the EM algorithm
 
EM <- function(y,K,theta_init){
  
  LL = c()
  K  <- length(theta_init$mu)
  mu <- theta_init$mu
  omega <- theta_init$omega
  diff <-  1
  iter <- 0
  while (diff>10^(-6)){
    iter <- iter +1 
    #----------- E step
    tau <- sapply(1:K,function(k){dpois(y,mu[k])*omega[k]})
    norm <- rowSums(tau) %*%matrix(1,nrow=1,ncol=K)
    tau <- tau/norm
    
    #--------- M step 
    mu_new <- c(y %*% tau)/colSums(tau)
    omega_new <- colMeans(tau)
    # --- Step crit
    diff = sum((mu_new - mu)^2) 
    mu <- mu_new 
    omega <- omega_new
    #---- LL
    LLiter <- LogLikPoissonMixture(y,mu,omega)
    LL <- c(LL,LLiter )
    print(c(iter,LLiter ))
          
  }
  return(list(LL = LL,mu=mu,omega=omega))
}
 


### 5  Initialisation 
init = kmeans(y,seq(1:K))
theta_init <- list(mu = c(init$centers), omega=init$size/n)

### 6 run EM 

res_EM <-EM(y,K,theta_init)
plot(res_EM$LL,type='l')




### 7 For several K 

chooseK <- function(y,Kmax){
 
  LL_K <- rep(0,Kmax) 
  LL_K[1] =sum( dpois(y,mean(y),log = TRUE))
  res_EM = list()
  res_EM[[1]] = list(mu=mean(y),omega=1,LL = LL_K[1],tau = rep(1,n))
  for (K in 2:Kmax){
  print(K)
  init = kmeans(y,seq(1:K))
  theta_init <- list(mu = c(init$centers), omega=init$size/n)
  res_EM[[K]] <-EM(y,K,theta_init)
  LL_K[K] <- max(res_EM[[K]]$LL)
  }
  NbParam <- 2*(1:Kmax)-1
  AIC = LL_K - NbParam
  BIC = LL_K - log(n)/2*NbParam
  return(list(LL_K=LL_K,AIC = AIC,BIC = BIC))
}

resSelectK <- chooseK(y,Kmax) 

plot(1:Kmax,resSelectK$LL_K,type='l')
lines(1:Kmax,resSelectK$BIC,col='red'); abline(v = which.max(resSelectK$BIC),col='red',lty=2)
lines(1:Kmax,resSelectK$AIC,col='green');abline(v = which.max(resSelectK$AIC),col='green',lty=3)



########## 8 # change the Data

Z = sample(1:3,n,replace= TRUE,prob   = omega)
lambdaZ = rgamma(n,0.4*mu[Z],0.4)
yPG = rpois(n,lambdaZ)
hist(yPG)

Kmax = 10
resSelectK <- chooseK(yPG,Kmax) 

plot(1:Kmax,resSelectK$LL_K,type='l')
lines(1:Kmax,resSelectK$BIC,col='red'); abline(v = which.max(resSelectK$BIC),col='red',lty=2)
lines(1:Kmax,resSelectK$AIC,col='green');abline(v = which.max(resSelectK$AIC),col='green',lty=3)


