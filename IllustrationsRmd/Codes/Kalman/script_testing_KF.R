library(tidyverse)
source("script_simulate_data.R")
source("utils_kalman_functions.R")
pars0 <- list(A = diag(1, 2), B = diag(0.8, 2), Sigma = diag(1, 2),
              Omega = diag(0.1, 2))
pars_true <- list(A = A_true, B = B_true, Sigma = Sigma_true,
              Omega = Omega_true)
filt_res  <- get_filter_quants(pars_true, Y)
smooth_res <- get_smoother_quants(pars_h = pars_true, 
                                  filt_pars = filt_res)

n_steps <- 10000
for(i in 1:n_steps){
  pars0$Sigma <- Sigma_true
  pars0$Omega <- Omega_true
  filt_res <- get_filter_quants(pars0, Y)
  smooth_res <- get_smoother_quants(pars_h = pars0, 
                                    filt_pars = filt_res)
  pars0 <- get_new_pars(smooth_res, Y)
  print(pars0$A %*% pars0$B)
}

