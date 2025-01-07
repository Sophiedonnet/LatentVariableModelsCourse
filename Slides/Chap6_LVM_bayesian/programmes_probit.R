rm(list=ls())
library(mvtnorm)
library(coda)
############### PROBIT


### simulation des données
n=240
X= matrix(1,n,3)
X[,2] = floor(rgamma(n,90*0.5,1*0.5))
X[,3] = floor(X[,2]*runif(n,1/2,1))
U = cov(X)
X[,2] = (X[,2]-mean(X[,2]))/U[2,2]
X[,3] = (X[,3]-mean(X[,3]))/U[3,3]
theta_vrai=c(0,4,-5)
p = pnorm(X%*%matrix(theta_vrai,ncol=1))
Y = rbinom(n,rep(1,n),p)
data_alz = list(Y=Y,age=X[,2],debut = X[,3])

###  GLM 
res_glm = glm(Y~age+debut, family = binomial(link = "probit"), data=data_alz)
theta_chap=res_glm$coefficients
Sigma = summary(res_glm)$cov.unscaled
#### 

### MH
# initialisation
theta=theta_chap

# paramètres 

MH_probit = function(Y,X,theta,B,tau,Sigma){ 
  
  ### Metropolis_hastings pour un modèle de régression probit. 
  
  Theta_stock = matrix(0,B,3)
  log_lik = sum(dbinom(Y,rep(1,n),pnorm(X%*%matrix(theta,ncol=1)),log=TRUE))
  LL=rep(0,B)
  
  acc=0
  for (b in 1:B){
    if(b%%5000==0){print(b)}
    theta_c = theta+rmvnorm(1,rep(0,3),tau^2*Sigma)
    log_lik_c = sum(dbinom(Y,rep(1,n),pnorm(X%*%matrix(theta_c,ncol=1)),log=TRUE))
    rho = log_lik_c-log_lik; 
    U=runif(1)
    if (log(U)<rho){
      acc=acc+1
      theta=theta_c
      log_lik = log_lik_c
    }
    Theta_stock[b,] =theta
    LL[b] = log_lik
  }
  res=list(LL=LL,Theta_stock=Theta_stock,acc=acc/B)
  return(res)
}

Gibbs_probit = function(Y,X,theta,B){
  p = length(theta)
  Theta_stock = matrix(0,B,p)
  LL=rep(0,B)
  invtXX = solve(t(X)%*%X)
  
  for (b in 1:B){
    if(b%%5000==0){print(b)}
    # Z | theta, Y  : gaussiennes tronquées
    mu = X%*%matrix(theta,ncol=1)
    U=runif(n)
    Z = qnorm(U*((Y==1)-(Y==0))*pnorm(mu)+U*(Y==0)+(1-pnorm(mu))*(Y==1))+mu
    # theta | Z, Y  : gaussiennes tronquées
    
    theta = rmvnorm(1,invtXX%*%t(X)%*%matrix(Z,ncol=1),invtXX)
    Theta_stock[b,] =theta
    LL[b] = sum(dbinom(Y,rep(1,n),pnorm(X%*%matrix(theta,ncol=1)),log=TRUE))
    
  }
  res=list(LL=LL,Theta_stock=Theta_stock)
  return(res)
}




B=50000

### 3 metropolis Hastings
tau_1 = 0.01
res_MH_1  = MH_probit(Y,X,theta,B,tau_1,Sigma)
Theta_1 = res_MH_1$Theta_stock
acc_1 = res_MH_1$acc

tau_2 = 10
res_MH_2  = MH_probit(Y,X,theta,B,tau_2,Sigma)
Theta_2 = res_MH_2$Theta_stock
acc_2 = res_MH_2$acc

tau_3 = 1.5
res_MH_3  = MH_probit(Y,X,theta,B,tau_3,Sigma)
Theta_3 = res_MH_3$Theta_stock
acc_3 = res_MH_3$acc

part=(B/2):B
THETA = list(Theta_1,Theta_2,Theta_3)
par(mfrow=c(3,3))
for (l in 1:3){
for(k in 1:3){
  plot(part,THETA[[l]][part,k],type='l',xlab="",ylab="",main="")
  abline(h=theta_vrai[k],col='red')
  }
  
}

acc_1
acc_2
acc_3
par(mfrow=c(3,3))

autocorr.plot(MH(Theta_1[part,1],start=25000))
autocorr.plot(MH(Theta_1[part,2],start=25000))
autocorr.plot(MH(Theta_1[part,3],start=25000))

autocorr.plot(MH(Theta_2[part,1],start=25000))
autocorr.plot(MH(Theta_2[part,2],start=25000))
autocorr.plot(MH(Theta_2[part,3],start=25000))

autocorr.plot(MH(Theta_3[part,1],start=25000))
autocorr.plot(MH(Theta_3[part,2],start=25000))
autocorr.plot(MH(Theta_3[part,3],start=25000))

################### 

res_Gibbs=Gibbs_probit(Y,X,theta,B)
par(mfrow=c(2,3))
for(k in 1:3){
  plot(part,res_MH_3$Theta_stock[part,k],type='l',xlab="",ylab="",main="")
  abline(h=theta_vrai[k],col='red')
}
for(k in 1:3){
  plot(part,res_Gibbs$Theta_stock[part,k],type='l',xlab="",ylab="",main="")
  abline(h=theta_vrai[k],col='red')
}

part2 = seq(B/2,B,by=6)
par(mfrow=c(1,3))
for(k in 1:3){
  plot(density(res_MH_3$Theta_stock[part2,k]),col='blue',lwd=2,main="",ylab="",xlab="")
  lines(density(res_Gibbs$Theta_stock[part2,k]),col='orange',lwd=2)
  abline(v=theta_vrai[k],col='red')}

autocorr.plot(mcmc(res_Gibbs$Theta_stock[part,1],start=25000))
autocorr.plot(mcmc(res_Gibbs$Theta_stock[part,2],start=25000))
autocorr.plot(mcmc(res_Gibbs$Theta_stock[part,3],start=25000))

