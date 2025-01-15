

n = 10;
p = 0.3;

X = 2#rbinom(1,n,p)

lik = function(p){
	p^{X}*(1-p)^(n-X)
}


curve(lik,0,1,col='red',lwd=2,xlab='p',ylab='Verosimilitud')
abline(v=X/n,lty=2)


pc = X/n;


A = pc - qnorm(0.975)*sqrt(pc*(1-pc)/n)

B = pc + qnorm(0.975)*sqrt(pc*(1-pc)/n)


absc =seq(-0.1,1.1,length=1000) 
Y.absc = dbeta(absc,X+1,n-X+1)
curve(dbeta(x,X+1,n-X+1),-0.1,1.1,col='red',lwd=2,xlab = 'Values of p', ylab ='Prior y posterior distributions') 

lines(absc,dunif(absc),col='blue',lwd=2) 


(X+1)/(n+2)

alpha = X+1
beta = n-X+1

b = qbeta(0.975,alpha,beta)
a = qbeta(0.025,alpha,beta)

alpha = X+1
beta = n-X+1




par(mfrow=c(1,1))
x  = seq(a,b,len=1000)
x.decrois = seq(b,a,len=1000)
y = dbeta(x.decrois,alpha,beta)

xx = c(x,rep(b,100),x.decrois,rep(a,100));
yy = c(rep(-0.5,1000),rep(dbeta(b,alpha,beta),100),y,rep(dbeta(a,alpha,beta),100))

p1 = alpha/(alpha+beta)
p2 = (alpha-1)/(alpha+beta-2)

par(cex = 1)
plot(absc, Y.absc, type = "l",lty=1,lwd=2,col='red', ylab = "Prior y posterior distributions", xlab = 'Values of p',main='')
lines(absc,dunif(absc),col='blue',lwd=2) 
polygon(xx, yy,density = c(5),col='red')
abline(v=a,col='red',lty=1,lwd=2)
abline(v=b,col='red',lty=1,lwd=2)
abline(v=p1,col='green',lwd=2)
abline(v=p2,col='black',lwd=2)
legend(0.7,3,col=c("red","green","black"),legend=c("High pobability interval","Posterior mode","Posterior mean"),lty = c(1,1,1),lwd=c(2,2,2))





absc =seq(-0.1,1.1,length=1000) 
abscneg = seq(-0.1,-0.001,length=10)
abscpos = seq(0,1.1,length=1000)
prior0_0 = dunif(absc)
prior1_3  = dbeta(absc[-(1:5)],1,3)
prior2_5  = dbeta(absc,2,5)

plot(absc,prior2_5,col='orange',type='l',lwd=2,ylab = "Prior distributions", xlab = 'Values of p',main='',ylim=c(0,3))
lines(abscpos,dbeta(abscpos,1,3),col='magenta',type='l',lwd=2)
lines(abscneg,dbeta(abscneg,1,3),col='magenta',lty=1,lwd=2)
lines(absc,prior0_0,col='blue',lty=1,lwd=2)

legend(0.7,3,col=c("blue","magenta","orange"),legend=c("Uniform prior","Beta(1,3)","Beta(2,5)"),lty = c(1,1,1),lwd=c(2,2,2))


absc =seq(-0.1,1.1,length=1000) 
abscneg = seq(-0.1,-0.001,length=10)
abscpos = seq(0,1.1,length=1000)
post0_0 = dbeta(absc,1+X,1+n-X)
post1_3  = dbeta(absc,1+X,3+n-X)
post2_5  = dbeta(absc,2+X,5+n-X)

plot(absc,post2_5,col='orange',type='l',lwd=2,ylab = "Posterior distributions", xlab = 'Values of p',main='',ylim=c(0,4.3))
lines(abscpos,dbeta(abscpos,1+X,3+n-X),col='magenta',type='l',lwd=2)
lines(abscneg,dbeta(abscneg,1+X,3+n-X),col='magenta',lty=1,lwd=2)
lines(absc,post0_0,col='blue',lty=1,lwd=2)
lines(absc,prior2_5,col='orange',lwd=1,lty=2)
lines(abscpos,dbeta(abscpos,1,3),col='magenta',lwd=1,lty=2)
lines(abscneg,dbeta(abscneg,1,3),col='magenta',lty=2,lwd=1)
lines(absc,prior0_0,col='blue',lty=2,lwd=1)



legend(0.7,3,col=c("blue","magenta","orange"),legend=c("Uniform prior","Prior Beta(1,3)","Prior Beta(2,5)"),lty = c(1,1,1),lwd=c(2,2,2))






 n1 = 10; X1 = 2
n2 = 100; X2  = rbinom(1,n2,p)
n3 = 1000; X3  = rbinom(1,n3,p)
n4 = 10000; X4  = rbinom(1,n4,p)

absc =seq(-0,0.6,length=1000) 
prior = dbeta(absc,2,5)
post1  = dbeta(absc,2+X1,5+n1-X1)
post2  = dbeta(absc,2+X2,5+n2-X2)
post3  = dbeta(absc,2+X3,5+n3-X3)
#post4  = dbeta(absc,2+X4,5+n4-X4)


par(mfrow = c(1,1))
plot(absc,post3,col='orange',type='l',lwd=2,ylab = "Posterior distributions", xlab = 'Values of p',main='',lty=2)
lines(absc,prior2_5,col='orange',lwd=1,lty=1)
lines(absc,post2,col='orange',lwd=2,lty=3)
#lines(absc,post3,col='orange',lty=2,lwd=1)
lines(absc,post1,col='orange',lty=4,lwd=2)



post1  = dbeta(absc,1+X1,3+n1-X1)
post2  = dbeta(absc,1+X2,3+n2-X2)
post3  = dbeta(absc,1+X3,3+n3-X3)
#plot(absc,post3,col='magenta',type='l',lty=2,lwd=1,ylab = "Posterior distributions", xlab = 'Values of p',main='')
lines(absc,post3,col="magenta",lty=2,lwd=2)
lines(abscpos,dbeta(abscpos,1,3),col='magenta',type='l',lwd=1)
lines(abscneg,dbeta(abscneg,1,3),col='magenta',lty=1,lwd=2)
lines(absc,post2,col='magenta',lwd=2,lty=3)
#lines(absc,post3,col='orange',lty=2,lwd=1)
lines(absc,post1,col='magenta',lty=4,lwd=2)


legend(0,20,col=rep(1,4),legend=c("Prior","n=10","n=100","n=1000"),lty = c(1,4,3,2),lwd=c(1,2,2,2))






