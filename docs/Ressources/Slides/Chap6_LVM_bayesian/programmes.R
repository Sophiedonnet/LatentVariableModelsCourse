curve(dbeta(x,1/2,1/2),0,1,ylab="Beta density",xlab="p",lwd=2,col=2,ylim=c(0,4))
curve(dbeta(x,1,1),0,1,lwd=2,col=1,add=TRUE)
curve(dbeta(x,2,2),0,1,lwd=2,col=3,add=TRUE)
curve(dbeta(x,2,4),0,1,lwd=2,col=4,add=TRUE)
curve(dbeta(x,6,4),0,1,lwd=2,col=5,add=TRUE)

legend(0.2,4, # places a legend at the appropriate place 
       c("Beta(1/2,1/2)","Beta(1,1) = uniforme","Beta(2,2)","Beta(2,4)","Beta(6,4)"), # puts text in the legend
       lty=c(1,1,1,1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2,2,2,2,2),col=c(2,1,3,4,5)) # gives the legend lines the correct color and width



##########################" 
a_k=1; b_k=1
library(BayesLCA)
data("Alzheimer")
Y=Alzheimer
n  = dim(Y)[1]
K = dim(Y)[2]
n_k = apply(Y,2,sum)
titre=c("Halluc.(p1)", "Tr. attention (p2)", "Agress.(p3)", "Agit.(p4)", "Perturb ryth. circ (p5)", "Tr. affectifs (p6)")
par(mfrow=c(2,3))
for (k in 1:K){
  curve(dbeta(x,1,1),0,1,ylim=c(0,20),lwd=2,ylab='Post',main=titre[k])
  curve(dbeta(x,n_k[k]+1,n-n_k[k]+1),col='red',lwd=2,add=TRUE)
}
####################################" 

##########################" 
par(mfrow=c(1,1))


n=0; n_1=0
curve(dbeta(x,1/2+n_1,1/2+n-n_1),0,1,ylab="",main='Prior',xlab="p",lwd=2,col=2,ylim=c(0,3))
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,col=1,add=TRUE)
curve(dbeta(x,2+n_1,2+n-n_1),0,1,lwd=2,col=3,add=TRUE)
curve(dbeta(x,2+n_1,4+n-n_1),0,1,lwd=2,col=4,add=TRUE)
curve(dbeta(x,6+n_1,4+n-n_1),0,1,lwd=2,col=5,add=TRUE)



n=10; n_1=sum(Y[1:n,1])
curve(dbeta(x,1/2+n_1,1/2+n-n_1),0,0.4,ylab="",xlab="p",lwd=2,col=2,ylim=c(0,20),main="n=10")
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,col=1,add=TRUE)
curve(dbeta(x,2+n_1,2+n-n_1),0,1,lwd=2,col=3,add=TRUE)
curve(dbeta(x,2+n_1,4+n-n_1),0,1,lwd=2,col=4,add=TRUE)
curve(dbeta(x,6+n_1,4+n-n_1),0,1,lwd=2,col=5,add=TRUE)


n=50; n_1=sum(Y[1:n,1])
curve(dbeta(x,1/2+n_1,1/2+n-n_1),0,0.4,ylab="",xlab="p",lwd=2,col=2,ylim=c(0,20),main="n=50")
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,col=1,add=TRUE)
curve(dbeta(x,2+n_1,2+n-n_1),0,1,lwd=2,col=3,add=TRUE)
curve(dbeta(x,2+n_1,4+n-n_1),0,1,lwd=2,col=4,add=TRUE)
curve(dbeta(x,6+n_1,4+n-n_1),0,1,lwd=2,col=5,add=TRUE)


n=240; n_1=sum(Y[1:n,1])
curve(dbeta(x,1/2+n_1,1/2+n-n_1),0,0.4,ylab="",xlab="p",lwd=2,col=2,ylim=c(0,20),main="n=240")
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,col=1,add=TRUE)
curve(dbeta(x,2+n_1,2+n-n_1),0,1,lwd=2,col=3,add=TRUE)
curve(dbeta(x,2+n_1,4+n-n_1),0,1,lwd=2,col=4,add=TRUE)
curve(dbeta(x,6+n_1,4+n-n_1),0,1,lwd=2,col=5,add=TRUE)


####################################" 



#---------------------------------------------------------------------------- 

library(gtools)

G=3; 
Znum=sample(1:G,n,rep=TRUE,prob=rep(1/G,G))
a.prior = 1;
b.prior = 1
d.prior = 1

pi = rep(1/G,G)
gamma = matrix((n_k+1)/(n_k+1+n-n_k+1),nrow=1)
B=10000; 

Y = as.matrix(Y)

Gibbs.MixtMultiBern <- function(Y, Znum, pi, gamma, a.prior, b.prior, d.prior, B = 5){
  
  # Y = data matrix: n * J
  # Znum = classification matrix : n * 1
  
  
  
  
  n = length(Znum); J = ncol(gamma); K = length(pi)
  Gamma_stock=array(0,c(B,K,J))
  Pi_stock=matrix(0,B,K)
  
  Z = vapply(1:K, function(k){as.numeric(Znum==k)}, rep(1,n))
  for (b in 1:B){
    if(b%%1000==0){print(b)}
    N = colSums(Z)
    pi = rdirichlet(1, N + d.prior)
    gamma = matrix(rbeta(K*J, t(Z)%*%Y + a.prior, t(Z)%*%(Y==0) + b.prior), K, J)    
    tau = exp(Y%*%t(log(gamma)) + (1-Y)%*%t(log(1-gamma)) + rep(1, n)%*%log(pi))
    tau = tau/rowSums(tau)
    Z = t(sapply(1:n, function(i){rmultinom(1, 1, tau[i, ])}))
    Gamma_stock[b,,]= gamma
    Pi_stock[b,]=pi
    
  }
  return(list(Z=Z, Znum = Z%*%(1:K), Gamma_stock=Gamma_stock, Pi_stock=Pi_stock))
}
RES = Gibbs.MixtMultiBern(Y, Znum, pi, gamma, a.prior, b.prior, d.prior, B=100000)



par(mfrow=c(G,1))
for (g in 1:G){
  plot(c(1,B), c(0,1),type='n')
  lines(RES$Pi_stock[,g],col=g)
}
  
  




############### PROBIT
n=240
X= matrix(1,n,3)
X[,2] = floor(rgamma(n,90*0.5,1*0.5))
X[,3] = floor(X[,2]*runif(n,1/2,1))

U = cov(X)
X[,2] = (X[,2]-mean(X[,2]))/U[2,2]
X[,3] = (X[,3]-mean(X[,3]))/U[3,3]


theta_vrai=c(0,4,-5)
q = X%*%matrix(theta_vrai,ncol=1); hist(q)
p = pnorm(X%*%matrix(theta_vrai,ncol=1))
hist(p)
Y = rbinom(n,rep(1,n),p)

data_alz = list(Y=Y,age=X[,2],debut = X[,3])
res_glm = glm(Y~age+debut, family = binomial(link = "probit"), data=data_alz)
theta_chap=res_glm$coefficients
Sigma = summary(res_glm)$cov.unscaled

library(mvtnorm)

theta=rep(-10,3)
B=50000

Theta_stock = matrix(0,B,3)
log_lik = sum(dbinom(Y,rep(1,n),pnorm(X%*%matrix(theta,ncol=1)),log=TRUE))
LL=rep(0,B)
tau=1.3
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
acc/B
par(mfrow=c(2,2))
part=(1):B
plot(LL,type='l')
for (i in 1:3){
  plot(Theta_stock[part,i],type='l')
  abline(h=theta_chap[i],col='red')
  abline(h=theta_vrai[i],col='green')
}
  
res_mcmc=mcmc(Theta_stock)
autocorr.plot(res_mcmc)
gelman.diag(res_mcmc)
