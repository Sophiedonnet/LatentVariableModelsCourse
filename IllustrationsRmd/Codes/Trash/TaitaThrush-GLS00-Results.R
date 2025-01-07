# Results of genotype mixture and admixture for Taita thrush data (GLS00)

rm(list=ls()); palette('R3'); par(pch=20)

# -> Book: TaitaThrush-k5-m8-5-6-3-4-10-8-seed1-try20.Rdata

# Dirs
library(extraDistr); library(tensor); library(MSCquartets); library(sirt); 
library(ade4); library(mclust)
source('../Functions/FunctionsGenotypeMixture.R')
source('../Functions/FunctionsLatex.R')
dataDir <- '../../Data/TaitaThrush/'
dataName <- 'TaitaThrush'
figDir <- '../../Figures/'
exportFig <- FALSE

# Data & fit
kMax <- 5
dataName <- 'TaitaThrush'
fitName <- '-mixture-k5-m8-5-6-3-4-10-8-seed1-try20'
load(paste0(dataDir, dataName, fitName, '.Rdata'))
n <- nrow(genotypes); p <- ncol(genotypes);

################################################################################
# Mixture
alleles <- sapply(1:p, function(j){sort(unique(haplotypes[, j]))})
m <- sapply(alleles, length)
logL <- sapply(1:kMax, function(k){emList[[k]]$logLik})
ent <- sapply(1:kMax, function(k){-sum(emList[[k]]$tau*log(emList[[k]]$tau))})
penDim <- (0:(kMax-1)) + (1:kMax)*(sum(m-1))
bic <- logL - penDim*log(n)/2
bicMax <- max(bic)
bicProb <- exp(bic-bicMax) / sum(exp(bic-bicMax))
icl <- bic - ent
rbind(logL, penDim, ent, bic, icl, bicProb)
if(exportFig){png(paste0(figDir, dataName, '-mixture-BIC.png'))}
plot(logL, type='b', ylim=range(c(logL, bic, icl)), pch=0, 
     ylab='log-likelihood', xlab=expression(K), cex=1.5, cex.axis=1.5, cex.lab=1.5)
lines(bic, type='b', col=2, pch=1, cex=1.5); abline(v=which.max(bic), col=2)
lines(icl, type='b', col=4, pch=2, cex=1.5); abline(v=which.max(icl), col=4)
if(exportFig){dev.off()}

# Fit
k <- which.max(bic)
em <- emList[[k]]
par(mfrow=c(ceiling(sqrt(p)), round(sqrt(p))), mex=.6)
dimnames(em$gamma)[[1]] <- 1:k
for(j in 1:p){
  if(exportFig){png(paste0(figDir, dataName, '-mixture-marker', j, '-gamma.png'))}
  barplot(t(em$gamma[, j, 1:m[j]]), main=colnames(genotypes)[j], 
          col=1+(1:m[j]),
          # col=gray.colors(m[j]), 
          cex=1.5, cex.axis=1.5)
  if(exportFig){dev.off()}
}

# Clustering
par(mfrow=c(1, 1), pch=20)
# bp <- barplot(t(em$pi[order(pop), ]), border=1+(1:k), col=0)
bp <- barplot(t(em$tau[order(pop), ]), col=4+(1:k), border=4+(1:k))
points(bp, rep(-.005, n), col=sort(pop), pch=14+sort(pop))
abline(v=.5*mean(diff(bp))+bp[cumsum(table(pop))][-4], col=1, lwd=3, lty=2)

if(k==3){
  simplexPrepare()
  for(i in 1:n){simplexPoint(em$tau[i, ], col=pop[i], pch=14+pop[i], cex=1.5)}
  simplexPoint(em$pi, cex=3, pch='+', col=8)
}
Zhat <- apply(em$tau, 1, which.max)
table(pop, Zhat)

################################################################################
# Fit admixture (with fixed prob_i)
kMax <- 5; tryNb <- 20
resName <- paste0(dataName, '-admixture-k', kMax, '-m')
for(j in 1:p){resName <- paste0(resName, m[j], '-')}
algoParms <- paste0('seed', seed, '-try', tryNb)
resFile <- paste0(dataDir, resName, algoParms, '.Rdata')
if(!file.exists(resFile)){
  emList <- initList <- list(); par(mfrow=c(ceiling(sqrt(2*kMax)), round(sqrt(2*kMax))))
  for(k in 1:kMax){
    # Init one tau / marker
    cat('\n k =', k, 'each:')
    tauInit <- array(dim=c(n, p, k))
    for(j in 1:p){tauInit[, j, ] <- InitEM(genotypes=genotypes[, j, drop=FALSE], k=k)$tau}
    if(k > 1){for(j in 1:p){
      tauInit[, j, ] <- tauInit[, j, ] + 1/k
      tauInit[, j, ] <- tauInit[, j, ] / rowSums(tauInit[, j, ])
    }}
    emList[[k]] <- EMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                       k=k, tauInit=tauInit)
    cat(' ->', emList[[k]]$logLik)
    plot(emList[[k]]$logLikPath, type='b', ylab='', xlab='', main=paste('k =', k), 
         ylim=quantile(emList[[k]]$logLikPath, prob=c(.1, 1)))
    # Init same tau
    cat('\n k =', k, 'same:')
    tauInit <- array(dim=c(n, p, k))
    tauInit[, 1, ] <- InitEM(genotypes=genotypes[, j, drop=FALSE], k=k)$tau
    for(j in 2:p){tauInit[, j, ] <- tauInit[, 1, ]}
    if(k > 1){for(j in 1:p){
      tauInit[, j, ] <- tauInit[, j, ] + 1/k
      tauInit[, j, ] <- tauInit[, j, ] / rowSums(tauInit[, j, ])
    }}
    emTmp <- EMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                         k=k, tauInit=tauInit)
    cat(' ->', emTmp$logLik)
    if(emTmp$logLik > emList[[k]]$logLik){cat(' !!!'); emList[[k]] <- emTmp}
    # Random init one tau / marker
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try each', tt, ':')
      tauTry <- array(dim=c(n, p, k))
      for(j in 1:p){tauTry[, j, ] <- rdirichlet(n, rep(1, k))}
      emTmp <- EMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                     k=k, tauInit=tauTry)
      cat(' ->', emTmp$logLik)
      if(emTmp$logLik > emList[[k]]$logLik){cat(' !!!'); emList[[k]] <- emTmp}
    }}
    # Random init same tau / marker
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try same', tt, ':')
      tauTry <- array(dim=c(n, p, k))
      tauTry[, 1, ] <- rdirichlet(n, rep(1, k))
      for(j in 2:p){tauTry[, j, ] <- tauTry[, 1, ]}
      emTmp <- EMgenotypeAdMixtureHW(haplotypes=haplotypes, haploArray=haploArray, 
                                     k=k, tauInit=tauTry)
      cat(' ->', emTmp$logLik)
      if(emTmp$logLik > emList[[k]]$logLik){cat(' !!!'); emList[[k]] <- emTmp}
    }}
    plot(emList[[k]]$logLikPath, type='b', ylab='', xlab='', main=paste('k =', k), 
         ylim=quantile(emList[[k]]$logLikPath, prob=c(.1, 1)))
  }
  save(haplotypes, haploArray, genotypes, initList, emList, file=resFile)
}else{load(resFile)}
  
# Dirichlet for the pi_i
dirList <- list()
dirList[[1]] <- list(alpha=1, alpha0=1, xsi=1, logLik=0)
for(k in 2:kMax){
  dirList[[k]] <- dirichlet.mle(emList[[k]]$pi)
  dirList[[k]]$logLik <- sum(sapply(1:n, function(i){ddirichlet(emList[[k]]$pi[i, ], dirList[[k]]$alpha, log=TRUE)}))
  }

# Results
par(mfrow=c(1, 1))
alleles <- sapply(1:p, function(j){sort(unique(haplotypes[, j]))})
m <- sapply(alleles, length)
logL <- sapply(1:kMax, function(k){emList[[k]]$logLik})
penDim <- (0:(kMax-1)) + (1:kMax)*(sum(m-1))
# logL <- sapply(1:kMax, function(k){emList[[k]]$logLik + dirList[[k]]$logLik})
# penDim <- (1:kMax) + (1:kMax)*(sum(m-1))
ent <- sapply(1:kMax, function(k){-sum(emList[[k]]$tau*log(emList[[k]]$tau))})
bic <- logL - penDim*log(n)/2
bicMax <- max(bic)
bicProb <- exp(bic-bicMax) / sum(exp(bic-bicMax))
icl <- bic - ent
rbind(logL, penDim, ent, bic, icl, bicProb)
plot(logL, type='b', ylim=range(c(logL, bic, icl)))
lines(bic, type='b', col=2); abline(v=which.max(bic), col=2)
lines(icl, type='b', col=4); abline(v=which.max(icl), col=4)

# Fit
k <- which.max(bic)
k <- 3
em <- emList[[k]]
par(mfrow=c(ceiling(sqrt(p)), round(sqrt(p))), mex=.6)
for(j in 1:p){
  barplot(t(em$gamma[, j, 1:m[j]]), main=colnames(genotypes)[j], col=1+(1:m[j]))
}

# Clustering
par(mfrow=c(1, 1), pch=20)
# bp <- barplot(t(em$pi[order(pop), ]), border=1+(1:k), col=0)
bp <- barplot(t(em$pi[order(pop), ]), col=4+(1:k), border=4+(1:k))
points(bp, rep(-.005, n), col=sort(pop), pch=14+sort(pop))
abline(v=.5*mean(diff(bp))+bp[cumsum(table(pop))][-4], col=1, lwd=3, lty=2)
simplexPrepare()
for(i in 1:n){simplexPoint(em$pi[i, ], col=pop[i], pch=14+pop[i], cex=1.5)}
simplexPoint(dirList[[k]]$xsi, cex=3, pch='+', col=8)
Zhat <- apply(em$pi, 1, which.max)
table(pop, Zhat)

