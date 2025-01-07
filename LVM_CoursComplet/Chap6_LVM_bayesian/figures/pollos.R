
A = exp(8)
B = 5 ; 
C = 14 /100
logf = function(t,log_A,log_B,log_C){
	log_A-exp(log_B)*exp(-exp(log_C)*t)}
	
	
	
t = c(0,2,4,6,8,12,16,20,24,28,32,36,40)
n = length(t)
sigma_vrai = 0.05
#y = logf(t,8,log(5),log(14/100))+rnorm(n,0,sigma_vrai)	
#plot(t,y)	
#curve(logf(x,8,log(5),log(14/100)),add=TRUE,col='red')	

#save(y,file="weight.R")
load("weight.R")





plot(t,exp(y),type='p',pch=18,xlab="Time (weeks)",ylab="Chick's weight")
lines(t,exp(y))


alpha = 0 ; 
beta =0;

mu_log_A = log(3000)
mu_log_B = 3
omega_log_A = 3
omega_log_B = 3


log_A = log(3000)
log_C = log(14/100)
log_B = log(5)

M = 400000
rho=0.1

sigma2 = rep(0,M)
log_A = rep(log(3000),M); log_A[1] =log(3000)
log_B = rep(log(5) ,M); log_B[1] = mu_log_B + 2
acc_B= 0
C = exp(log_C)


for (l in 2:M){
	
	if (l%%10000==0){	print(l)}
	###### simulation de sigma
	alpha_prime = n/2+alpha
	beta_prime = beta + 1/2*sum((y-logf(t,log_A[l-1],log_B[l-1],log_C))^2)
	sigma2[l] = 1/rgamma(1,shape=alpha_prime,rate=beta_prime)
	
	####### simulation de log A
	B_l = exp(log_B[l-1])
	
	omega_post_A = 1/(1/omega_log_A+n/sigma2[l])
	mean_post_A = omega_post_A*(mu_log_A/omega_log_A + 1/sigma2[l]*sum(y+B_l*exp(-C*t)));
	log_A[l] = rnorm(1,mean_post_A,sqrt(omega_post_A))
	
	####### simulation de log B
	v = sample(c(1,1/100,1/10))
	log_B_c = log_B[l-1]+rnorm(1,0,rho*v)
	
	p = sum(dnorm(y,logf(t,log_A[l],log_B_c,log_C),sqrt(sigma2[l])),log=TRUE) + dnorm(log_B_c,mu_log_B,omega_log_B,log=TRUE)-sum(dnorm(y,logf(t,log_A[l],log_B[l-1],log_C),sqrt(sigma2[l])),log=TRUE) - dnorm(log_B[l-1],mu_log_B,omega_log_B,log=TRUE)
	
	log_B[l] = log_B[l-1]
	if(log(runif(1))<p){
		
		log_B[l] = log_B_c
		acc_B = acc_B+1
	
	}
	}
taux_accept_B = acc_B/M
print(taux_accept_B*100)

par(mfrow=c(2,3))
plot(part,log_A[part],type='l')
abline(h = log(3000),col='red')
part = seq(2*M/4,M,by=1)
plot(part,log_B[part],type='l')
abline(h=log(5),col='red')
plot(part,sigma2[part],type='l')
abline(h = sigma_vrai^2,col='red')
acf(log_A[part],lag.max=100)
acf(log_B[part],lag.max=100)
acf(sigma2[part],lag.max=100)


quartz()
par(mfrow=c(1,3))
plot(density(log_A[part2]),main="log A",xlab="log_A")
curve(dnorm(x,mu_log_A,sqrt(omega_log_A)),add=TRUE,col='blue')
abline(v = log(3000),col='red')
plot(density(log_B[part2]),main="log B",xlab="log_B")
curve(dnorm(x,mu_log_B,sqrt(omega_log_B)),add=TRUE,col='blue')
abline(v = log(5),col='red')
plot(density(1/sigma2[part2]),main="1/sigma2",xlab="1/sigma2")
curve(dgamma(x,alpha,rate=beta),add=TRUE,col='blue')
abline(v=1/sigma_vrai^2,col='red')



quartz()
part2 = seq(M-50000,M,by=10);
reg = rep(0,length(t))
for (i in part2){
	reg = reg+logf(t,log_A[i],log_B[i],log_C)
}
reg = reg/length(part2)

plot(t,exp(y))
lines(t, exp(reg))
lines(t,exp(logf(t,8,log(5),log(14/100))),col='red',lty=2)



