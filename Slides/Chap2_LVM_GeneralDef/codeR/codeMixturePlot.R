mu = c(1,3,5)
sigma = c(1,1,3)
n = 1000
prop <- c(1/3,1/6,1/2)
Z = sample(c(1:3),n,replace=TRUE,prob = prop)
Y <- mu[Z] + sigma[Z]*rnorm(n)

hist(Y,freq = FALSE,nclass=100, main = '')
for (k in 1:3){
  curve(dnorm(x,mu[k],sigma[k])*prop[k],add=TRUE,col=k+1,lwd = 3)
}
abs = seq(-10,15,len=1000)
y <- dnorm(abs,mu[1],sigma[1])*prop[1]
for (k in 2:3){
  y <- y + dnorm(abs,mu[k],sigma[k])*prop[k]
}
lines(abs,y,col='black',lwd=3)