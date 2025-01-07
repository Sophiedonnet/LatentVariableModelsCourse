LogPhi <- function(genoTab, logGamma){sum(genoTab*logGamma, na.rm=TRUE)}

EMgenotypeMixture <- function(genotypes, genoArray, k, tauInit, tol=1e-6, iterMax=1e3, tauEps=1e-4){
  # Genotype mixture under 
  n <- nrow(genotypes); p <- ncol(genotypes)
  m <- sapply(1:p, function(j){length(unique(genotypes[, j]))})
  mMax <- max(m)
  tau <- tauInit
  diff <- 2*tol; iter <- 0
  logLik <- rep(NA, iterMax)
  while(diff > tol){
    iter <- iter+1; if(iter%%round(sqrt(iterMax))==0){cat('', iter)}
    # M step
    nGroup <- colSums(tau)
    pi <- nGroup/n
    logPi <- log(pi)
    gamma <- array(dim=c(k, p, mMax))
    for(g in 1:k){for(j in 1:p){for(a in 1:m[j]){
      gamma[g, j, a] <- sum(genoArray[, j, a]*tau[, g]/nGroup[g])
    }}}
    logGamma <- log(gamma)
    logPhi <- matrix(0, n, k)
    for(i in 1:n){for(g in 1:k){logPhi[i, g] <- LogPhi(genoArray[i, , ], logGamma[g, , ])}}
    
    # E step
    logTau <- logPhi
    logTau <- logTau + rep(1, n)%o%logPi
    logTau <- logTau - apply(logTau, 1, max)
    tauNew <- exp(logTau); tauNew <- tauNew / rowSums(tauNew)
    tauNew <- tauNew + tauEps; tauNew <- tauNew / rowSums(tauNew)
    
    # Test
    diff <- max(abs(tauNew - tau))
    tau <- tauNew
    logLik[iter] <- sum(tau%*%logPi) + sum(tau*logPhi) - sum(tau*log(tau))
  }
  cat('\n')
  return(list(tau=tau, logPhi=logPhi, gamma=gamma, pi=pi, m=m, 
              iter=iter, logLik=logLik[iter], logLikPath=logLik[1:iter]))
}
