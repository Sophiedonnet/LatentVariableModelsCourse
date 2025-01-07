# Test of MLE estimation for the Dirichlet multinomial model

rm(list=ls())
source('../Functions/FunctionsGenotypeMixture.R')
library(extraDistr); library(dirmult)

# Parms
n <- 1e2; p <- 7; k <- 3
alpha <- 1:k

# Data
w <- rdirichlet(n, alpha)
x <- t(sapply(1:n, function(i){rmultinom(1, p, prob=w[i, ])}))

# MLE Dirichlet
fit <- dirichlet.mle(x=w)
mle <- DirMLE(alpha=alpha, x=w)
rbind(alpha, fit$alpha, mle$alpha)

# MLE DirMult
fit <- dirmult(data=x, init=alpha)
mle <- DirMultMLE(alpha=alpha, x=x)
rbind(alpha, fit$gamma, mle$alpha)
