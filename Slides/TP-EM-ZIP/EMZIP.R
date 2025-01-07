 # ZIP model for abundance data

rm(list=ls()); 

# Data
abundance <- read.table('BarentsFish.csv', sep=';', header=TRUE)
j <- 10; n <- nrow(abundance)
y <- abundance$Se_me #[, 4+j]
y0 <- 1*(y==0)
hist(y, breaks=100)

# Functions
LogLik <- function(prob_pi, mu, y){sum(log(prob_pi*(y==0) + (1-prob_pi)*dpois(y, mu)))}

# Init
prob_pi <- mean(y==0)
mu <- mean(y[y>0])

LogLik(prob_pi,mu,y) 
# EM for ZIP model
tol <- 1e-6; diff <- 2*tol; 
iterMax <- 1e3; iter <- 1
logL <- rep(NA, iterMax)
logL[iter] <- LogLik(pi, mu, y)
while((diff > tol) & (iter < iterMax)){
  iter <- iter+1
  # E step
  tau <- y0*pi/(pi + (1-pi)*exp(-mu))
  tau <- cbind(tau, (1-tau)); 
  # M step
  piNew <- mean(tau[, 1])
  muNew <- sum(tau[, 2]*y)/sum(tau[, 2])
  # Test & update
  diff <- max(abs(c(piNew, muNew)-c(pi, mu)))
  pi <- piNew; 
  mu <- muNew
  logL[iter] <- LogLik(pi, mu, y)
  cat(pi, mu, diff, logL[iter], '\n')
}
plot(logL[1:iter], type='b', pch=20)

cat('ZIP: pi=', pi, ', mu=', mu, ', logL=', logL[iter], '\n')
cat('Poisson: mu=', mean(y), ', logL=', sum(dpois(y, mean(y), log=TRUE)), '\n')

