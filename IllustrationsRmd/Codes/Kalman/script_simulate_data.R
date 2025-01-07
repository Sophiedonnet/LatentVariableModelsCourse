require(tidyverse)
set.seed(123)
theta <- pi /12
B_true <- 0.99 * matrix(c(cos(theta), sin(theta), - sin(theta), cos(theta)),
                   nrow = 2)
Omega_true <- rWishart(1, 3, diag(0.05, 2))[,,1]
get_zs <- function(n_step, B_){
  dz <- ncol(B_)
  out <- matrix(nrow = n_step + 1,
                ncol = dz)
  out[1, ] <- 0
  for(i in 1:n_step){
    out[i + 1, ] <- B_ %*% out[i, ] + mixtools::rmvnorm(1, rep(0, 2), Omega_true)[1, ]
  }
  return(out)
}

n  <- 1e3 
Z_true <- get_zs(n, B_true)
A_true <- matrix(c(1, -0.2, 0, 0.8),
            nrow = 2)
Sigma_true <- round(rWishart(1, 3, diag(0.5, 2))[,,1], 2)

Y <- t(A_true %*% t(Z_true)) + mixtools::rmvnorm(n + 1, rep(0, 2), Sigma_true)

rm(theta, get_zs, n)