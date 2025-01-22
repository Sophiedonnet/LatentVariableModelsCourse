
InitEM <- function(genotypes, k, tauEps=1e-4){
  n <- nrow(genotypes); p <- ncol(genotypes)
  genoFrame <- list()
  for(j in 1:p){genoFrame[[j]] <- as.factor(genotypes[, j])}
  mca <- dudi.acm(genoFrame, scannf=FALSE, nf=10)
  cluster <- kmeans(mca$li, centers=k, nstart=100, iter.max=100)$cluster
  tau <- matrix(tauEps, n, k)
  for(i in 1:n){tau[i, cluster[i]] <- tau[i, cluster[i]] + 1}
  tau <- tau / rowSums(tau)
  return(list(mca=mca, cluster=cluster, tau=tau))
}

LogPhiHW <- function(haploTab, logGammaHaplo){
  sum(haploTab[1, , ]*logGammaHaplo, na.rm=TRUE) +
    sum(haploTab[2, , ]*logGammaHaplo, na.rm=TRUE) +
    log(2)*sum(haploTab[1, , ]*haploTab[2, , ], na.rm=TRUE)
  }




EMgenotypeMixtureHW <- function(haplotypes, haploArray, k=ncol(tauInit), tauInit, tol=1e-6, iterMax=1e3, tauEps=1e-4){
  # Genotype mixture under Hardy-Weinberg equilibrium
  # k=3; tauInit=initList[[k]]$tau; tol=1e-6; iterMax=1e1; tauEps=1e-4
  n <- nrow(haplotypes)/2; p <- ncol(haplotypes)
  geno2haplo <- matrix(0, n, 2*n)
  for(i in 1:n){geno2haplo[i, 2*(i-1)+1] <- geno2haplo[i, 2*(i-1)+2] <- 1}
  m <- sapply(1:p, function(j){length(unique(haplotypes[, j]))})
  mMax <- max(m)
  tauGeno <- tauInit; 
  diff <- 2*tol; iter <- 0
  logLik <- rep(NA, iterMax)
  while((diff > tol) & (iter < iterMax)){
    iter <- iter+1; if(iter%%round(sqrt(iterMax))==0){cat('', iter)}
    # M step
    nGroup <- colSums(tauGeno)
    pi <- nGroup/n
    logPi <- log(pi)
    tauHaplo <- t(geno2haplo) %*% tauGeno  
    gamma <- array(dim=c(k, p, mMax))
    for(g in 1:k){for(j in 1:p){for(a in 1:m[j]){
      gamma[g, j, a] <- sum(haploArray[, j, a]*tauHaplo[, g])/(sum((!is.na(haplotypes[, j]))*tauHaplo[, g]))
    }}}
    logGamma <- log(gamma)
    logPhi <- matrix(0, n, k)
    for(i in 1:n){for(g in 1:k){
      logPhi[i, g] <- LogPhiHW(haploArray[2*(i-1)+(1:2), , ], logGamma[g, , ])}}
    
    # E step
    logTau <- logPhi + rep(1, n)%o%logPi
    logTau <- logTau - apply(logTau, 1, max)
    tauNew <- exp(logTau); tauNew <- tauNew / rowSums(tauNew)
    tauNew <- tauNew + tauEps; tauNew <- tauNew / rowSums(tauNew)
    
    # Test
    diff <- max(abs(tauNew - tauGeno))
    tauGeno <- tauNew
    logLik[iter] <- sum(tauGeno%*%logPi) + sum(tauGeno*logPhi) - sum(tauGeno*log(tauGeno))
    # cat(logLik[iter], diff, '\n')
  }
  return(list(tau=tauGeno, logPhi=logPhi, gamma=gamma, pi=pi, m=m, 
              iter=iter, logLik=logLik[iter], logLikPath=logLik[1:iter]))
}

###########################################################""
#--------------------------------------------------------
###########################################################
LogDirichlet <-function(alpha){sum(lgamma(alpha))-lgamma(sum(alpha))} 

EntropyDirichlet <-function(alpha){LogDirichlet(alpha) + 
    (sum(alpha)-length(alpha))*digamma(sum(alpha)) -
    sum((alpha-1)*digamma(alpha))} 

Geno2Genotype <- function(geno){
  if(sum(is.na(geno)) > 0){geno <- geno[-which(rowSums(is.na(geno))>0), ]}
  p2 <- ncol(geno); genotypes <- matrix('', nrow(geno), p2/2)
  miss <- geno <- matrix(as.numeric(geno), nrow=nrow(geno), ncol=p2)
  for(i in 1:nrow(geno)){for(j in 1:(p2/2)){
    tmp <- sort(c(geno[i, 2*(j-1)+1], geno[i, 2*(j-1)+2]))
    if(tmp[1]==0){miss[i, j] <- 1}else{miss[i, j] <- 0}
    genotypes[i, j] <- paste0(tmp[1], tmp[2]) 
  }}
  return(list(genotypes=genotypes, miss=miss))
}

MakeGenoArray <- function(markers){
  n <- nrow(markers); p <- ncol(markers)
  alleles <- lapply(1:p, function(j){sort(unique(as.factor(markers[, j])))})
  m <- sapply(1:p, function(j){length(unique(as.factor(markers[, j])))})
  mMax <- max(m)
  markArray <- array(0, dim=c(n, p, mMax))
  for(i in 1:n){for(j in 1:p){
    markArray[i, j, which(alleles[[j]]==markers[i, j])] <- 1
  }}
  return(markArray)
}



