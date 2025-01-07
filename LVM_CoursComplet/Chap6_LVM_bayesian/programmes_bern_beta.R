### prior uniforme

par(mfrow=c(2,2))

Y=rbinom(500,1,0.4)

curve(dunif(x,0,1),0,1,ylab="",main='n=10',xlab="p",lwd=2,col=2,ylim=c(0,3))
n=10; n_1=sum(Y[1:n])
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,col=1,add=TRUE)

n=50; n_1=sum(Y[1:n])
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,lty=1,col=1,main='n=50',xlab="p",ylab="")
curve(dunif(x,0,1),0,1,ylab="",col=2,add=TRUE,lwd=2)



n=250; n_1=sum(Y[1:n])
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,lty=1,col=1,main='n=250',xlab="p",ylab="")
curve(dunif(x,0,1),0,1,ylab="",col=2,add=TRUE,lwd=2)



n=500; n_1=sum(Y[1:n])
curve(dbeta(x,1+n_1,1+n-n_1),0,1,lwd=2,lty=1,col=1,main='n=500',xlab="p",ylab="")
curve(dunif(x,0,1),0,1,ylab="",col=2,add=TRUE,lwd=2)

















library(BayesLCA)





data("Alzheimer")


curve(dbeta(x,1/2,1/2),0,1,ylab="Beta density",xlab="p",lwd=2,col=2,ylim=c(0,4))
curve(dbeta(x,1,1),0,1,lwd=2,col=1,add=TRUE)
curve(dbeta(x,2,2),0,1,lwd=2,col=3,add=TRUE)
curve(dbeta(x,2,4),0,1,lwd=2,col=4,add=TRUE)
curve(dbeta(x,6,4),0,1,lwd=2,col=5,add=TRUE)

legend(0.2,4, # places a legend at the appropriate place 
       c("Beta(1/2,1/2)","Beta(1,1) = uniforme","Beta(2,2)","Beta(2,4)","Beta(6,4)"), # puts text in the legend
       lty=c(1,1,1,1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2,2,2,2,2),col=c(2,1,3,4,5)) # gives the legend lines the correct color and width


Y=Alzheimer
n  = dim(Y)[1]
K = dim(Y)[2]
n_k = apply(Y,2,sum)

# ##########################" 


# a_k=1; b_k=1
# titre=c("Halluc.(p1)", "Tr. attention (p2)", "Agress.(p3)", "Agit.(p4)", "Perturb ryth. circ (p5)", "Tr. affectifs (p6)")
# par(mfrow=c(2,3))
# for (k in 1:K){
#   curve(dbeta(x,1,1),0,1,ylim=c(0,20),lwd=2,ylab='Post',main=titre[k])
#   curve(dbeta(x,n_k[k]+1,n-n_k[k]+1),col='red',lwd=2,add=TRUE)
# }
# ####################################" 



##########################" 



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






###################### INTERVALLES de Crédibilité ##" 
a=1
b=1
library(BayesLCA)
data("Alzheimer")
Y=Alzheimer
n  = dim(Y)[1]
K = dim(Y)[2]
n_k = apply(Y,2,sum)
n_1 = n_k[1]

a_post = a+n_1+1
b_post= b+n-n_1+1

################################################### 



q1 = qbeta(0.975,a_post,b_post)
q2 = qbeta(0.025,a_post,b_post)
absc =seq(0,0.2,length=1000) 
d.absc = dbeta(absc,a_post,b_post)



frame()






p1 = a_post/(a_post+b_post)
p2 = (a_post-1)/(a_post+b_post-2)



par(cex = 1)
k=1; 
plot(absc, d.absc, type = "l",lty=2,lwd=2,col='black', ylab="", xlab = 'Values of theta',main="")
q1 = 0.0451; q2 = 0.135

abline(h=k,col='red',lwd=2)
abline(v=q1,col='red',lty=2)
abline(v=q2,col='red',lty=2)

x  = seq(q1,q2,len=1000)
x.decrois = seq(q2,q1,len=1000)
y = dbeta(x.decrois,a_post,b_post)
xx = c(x,rep(q2,100),x.decrois,rep(q1,100));
yy = c(rep(-0.5,1000),rep(dbeta(q2,a_post,b_post),100),y,rep(dbeta(q1,a_post,b_post),100))
polygon(xx, yy,density = c(10),angle=45,lty=1,col='red',lwd=2)

prob1=pbeta(q2,a_post,b_post)-pbeta(q1,a_post,b_post)


par(cex = 1)
k=2;
plot(absc, d.absc, type = "l",lty=2,lwd=2,col='black', ylab = "", xlab = 'Values of theta',main=paste('k=',k,sep=''))
 
q1 = 0.0495; q2 = 0.127
abline(h=k,col='blue',lwd=2)
abline(v=q1,col='blue',lty=2)
abline(v=q2,col='blue',lty=2)

x  = seq(q1,q2,len=1000)
x.decrois = seq(q2,q1,len=1000)
y = dbeta(x.decrois,a_post,b_post)
xx = c(x,rep(q2,100),x.decrois,rep(q1,100));
yy = c(rep(-0.5,1000),rep(dbeta(q2,a_post,b_post),100),y,rep(dbeta(q1,a_post,b_post),100))
polygon(xx, yy,density = c(10),col='blue',lwd=2)

prob2=pbeta(q2,a_post,b_post)-pbeta(q1,a_post,b_post)



par(cex = 1)
k=3
plot(absc, d.absc, type = "l",lty=2,lwd=2,col='black', ylab = "", xlab = 'Values of theta',main='')

abline(h=k,col='green',lwd=2)
q1 = 0.0511; q2 = 0.123
prob3=pbeta(q2,a_post,b_post)-pbeta(q1,a_post,b_post)
abline(v=q1,col='green',lty=2,lwd=2)
abline(v=q2,col='green',lty=2,lwd=2)
x  = seq(q1,q2,len=1000)
x.decrois = seq(q2,q1,len=1000)
y = dbeta(x.decrois,a_post,b_post)
xx = c(x,rep(q2,100),x.decrois,rep(q1,100));
yy = c(rep(-0.5,1000),rep(dbeta(q2,a_post,b_post),100),y,rep(dbeta(q1,a_post,b_post),100))
polygon(xx, yy,density = c(10),angle=45,lty=1,col='green',lwd=2)


q1_et = qbeta(0.5*(1-0.9576843),a_post,b_post)
q2_et = qbeta(1-0.5*(1-0.9576843),a_post,b_post)
abline(v=q1_et,col='orange',lty=2,lwd=2)
abline(v=q2_et,col='orange',lty=2,lwd=2)



round(c(prob1*100,prob2*100,prob3*100),3)
f=function(x){dbeta(x,a_post/10,b_post/10)*1/2+dbeta(x,b_post/10,a_post/10)*1/2}
frame()
par(mfrow=c(1,1))
absc=seq(0,1,len=1000)
d.absc=f(absc)
plot(absc, d.absc, type = "l",lty=2,lwd=2,col='black', ylab = "", xlab = 'Values of theta',main="")
abline(h=1,col='red',lwd=2)
q1 = 0.007; q2 = 0.168;q3 = 0.833; q4=0.995;
abline(v=c(q1,q2,q3,q4),col='red',lty=2)
x  = seq(q1,q2,len=1000)
x.decrois = seq(q2,q1,len=1000)
y = f(x.decrois)
xx = c(x,rep(q2,100),x.decrois,rep(q1,100));
yy = c(rep(-0.5,1000),rep(f(q2),100),y,rep(f(q1),100))
polygon(xx, yy,density = c(10),angle=45,lty=1,col='red',lwd=2)

x  = seq(q3,q4,len=1000)
x.decrois = seq(q4,q3,len=1000)
y = f(x.decrois)
xx = c(x,rep(q4,100),x.decrois,rep(q3,100));
yy = c(rep(-0.5,1000),rep(f(q4),100),y,rep(f(q3),100))
polygon(xx, yy,density = c(10),angle=45,lty=1,col='red',lwd=2)




log(pbeta(0.1,a_post,b_post)/pbeta(0.1,1,1)/(1-pbeta(0.1,a_post,b_post))*(1-pbeta(0.1,1,1)),10)


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

  

