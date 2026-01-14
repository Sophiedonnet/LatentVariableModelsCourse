# ZIP model for abundance data

rm(list=ls()); 

# Data
abundance <- read.table('BarentsFish.csv', sep=';', header=TRUE)
j <- 10; n <- nrow(abundance)
y <- abundance$Se_ma
#y <- abundance$Tr_es
y0 <- 1*(y==0)
hist(y, breaks=sqrt(n))


# 
lik = function(y,mu,prob.pi){(y==0)*(1-prob.pi) + (prob.pi)*dpois(y,mu)}

#
prob.pi <- 1-mean(y==0)
mu <- mean(y[y>0])
log.lik <-c()
tol <- 1e-6; diff <- 2*tol; 
# 
while (diff > tol) {
  
  ### E step
  tau0 <- (1-prob.pi)*(y==0)/lik(y,mu,prob.pi)
  tau1 <- 1-tau0
  
  #### M step
  mu.new <- sum(y*tau1)/sum(tau1)
  prob.pi.new <- mean(tau0)
  
  #### diff mu, prob.pi
  diff = max(abs(mu.new-mu), abs(prob.pi.new-prob.pi)) 
  print(diff)
  mu <- mu.new
  prob.pi <- prob.pi.new
  
  ll<- sum(log(lik(y,mu,prob.pi)))
  print((c(mu,prob.pi,ll)))
  #####  log.Lik
  log.lik = c(log.lik,ll)
  
}
plot(log.lik,type='l')
