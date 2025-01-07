get_filter_quants <- function(pars_h, Y){
  # Model infos
  dz <- ncol(pars_h$B)
  dy <- ncol(Y)
  n_obs <- nrow(Y)
  V0 <- diag(1, dz)
  mu0 <- rep(0, dz)
  Iz <- diag(1, dz)
  # Creation of output objects
  filt_means <- matrix(nrow = n_obs, 
                       ncol = dz)
  filt_vars <- array(dim = c(dz, dz, n_obs))
  pred_vars <- array(dim = c(dz, dz, n_obs)) # W_t
  # Parameter extraction
  A_h <- pars_h$A; B_h <- pars_h$B; 
  Sigma_h <- pars_h$Sigma; Omega_h <- pars_h$Omega
  # Transformed parameters
  A_h_t <- t(A_h)
  B_h_t <- t(B_h)
  # Starting the filter
  ## Initialization
  Kgain <- V0 %*% A_h_t %*% solve((Sigma_h + A_h %*% V0 %*% A_h_t))
  filt_means[1, ] <- mu0 + Kgain %*% (Y[1, ] - A_h %*% mu0)
  filt_vars[,, 1] <- (Iz - Kgain %*% A_h) %*% V0
  for(i in 2:n_obs){
    pred_vars[,, i] <- Omega_h + B_h %*% filt_vars[,, i-1] %*% B_h_t
    Kgain <- pred_vars[,, i] %*% A_h_t %*% 
      solve(Sigma_h + A_h %*% pred_vars[,, i] %*% A_h_t)
    filt_means[i, ] <- B_h %*% filt_means[i-1, ] + 
      Kgain %*% (Y[i, ] - A_h %*% B_h %*% filt_means[i-1, ])
    filt_vars[,, i] <- (Iz - Kgain %*% A_h) %*% pred_vars[,, i]
  }
  list(filt_means = filt_means,
       filt_vars = filt_vars,
       pred_vars = pred_vars)
}



get_smoother_quants <- function(pars_h, filt_pars){
  # Model infos
  dz <- ncol(pars_h$B)
  dy <- ncol(Y)
  n_obs <- nrow(Y)
  V0 <- diag(1, dz)
  mu0 <- rep(0, dz)
  Iz <- diag(1, dz)
  # Creation of output objects
  smooth_means <- matrix(nrow = n_obs, 
                         ncol = dz)
  smooth_vars <- array(dim = c(dz, dz, n_obs))
  smooth_Ss <- array(dim = c(dz, dz, n_obs))
  smooth_covs <- array(dim = c(dz, dz, n_obs - 1)) # W_t
  smooth_Cs <- array(dim = c(dz, dz, n_obs - 1)) # W_t
  # Parameter extraction
  filt_means <- filt_pars$filt_means
  filt_vars <- filt_pars$filt_vars
  pred_vars <- filt_pars$pred_vars # W_t
  A_h <- pars_h$A; B_h <- pars_h$B; 
  Sigma_h <- pars_h$Sigma; Omega_h <- pars_h$Omega
  # Transformed parameters
  A_h_t <- t(A_h)
  B_h_t <- t(B_h)
  # Starting the filter
  ## Initialization
  smooth_means[n_obs, ] <- filt_means[n_obs, ]
  smooth_vars[,, n_obs] <- filt_vars[,, n_obs]
  smooth_Ss[,, n_obs] <- smooth_vars[,, n_obs] + 
    smooth_means[n_obs, ] %*% t(smooth_means[n_obs, ])
  for(i in (n_obs - 1):1){
    Ktilde <- filt_vars[,, i] %*% B_h_t %*% solve(pred_vars[,, i + 1])
    smooth_means[i, ] <- filt_means[i, ] + 
      Ktilde %*% (smooth_means[i + 1, ] - B_h %*% filt_means[i, ])
    smooth_vars[,, i] <- filt_vars[,, i] + 
      Ktilde %*% (smooth_vars[,, i+1] - pred_vars[,, i+1]) %*% t(Ktilde)
    smooth_covs[,, i] <- Ktilde %*% smooth_vars[,, i+1]
    smooth_Ss[,, i] <- smooth_vars[,, i] + 
      smooth_means[i, ] %*% t(smooth_means[i, ])
    smooth_Cs[,, i] <- smooth_covs[,, i] + 
      smooth_means[i, ] %*% t(smooth_means[i + 1, ])
  }
  list(smooth_means = smooth_means,
       smooth_vars = smooth_vars,
       smooth_covs = smooth_covs,
       smooth_Cs = smooth_Cs,
       smooth_Ss = smooth_Ss)
}



get_new_pars <- function(smooth_pars, Y){
  sum_Y_Yt = t(Y) %*% Y
  n_obs <- nrow(Y)
  sum_Y_mu <- t(Y) %*% smooth_pars$smooth_means
  sum_mu_Y <- t(sum_Y_mu)
  sum_S <- apply(smooth_pars$smooth_Ss,
                 c(1, 2), sum)
  Sum_S_m_n <- sum_S - smooth_pars$smooth_Ss[,,n_obs]
  Sum_S_m_1 <- sum_S - smooth_pars$smooth_Ss[,,1]
  sum_C <- apply(smooth_pars$smooth_Cs,
                 c(1, 2), sum)
  A <-   sum_Y_mu %*% solve(sum_S)
  # print("Sum_S_inv")
  # print(solve(sum_S))
  # print("Sum_Y_mu")
  # print(t(sum_Y_mu) )
  Sigma <- (sum_Y_Yt + A %*% sum_S %*% t(A) 
            - sum_Y_mu %*% t(A) - A %*% sum_mu_Y) / nrow(Y)
  B <- t(sum_C) %*%  solve(Sum_S_m_n)
  Omega <- (B %*%  Sum_S_m_n %*%  t(B) + Sum_S_m_1 
            - t(sum_C) %*% t(B) - B %*% sum_C) / (n_obs - 1)
  list(A = A, B = B, Sigma = Sigma, Omega = Omega)
}