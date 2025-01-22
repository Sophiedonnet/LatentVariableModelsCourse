# Fit genotype mixture & admixuture for Taita thrush data (GLS00)

rm(list=ls()); palette('R3'); par(pch=20)

# -> Book: TaitaThrush-k5-m8-5-6-3-4-10-8-seed1-try20.Rdata

# Dirs
 
library(ade4); 
library(extraDistr)
source('FunctionsGenotypeMixture.R')
dataDir <- '../../Data/TaitaThrush/'
dataName <- 'TaitaThrush'
 

# Data 
load('TaitaThrush-m8-5-6-3-4-10-8.Rdata')
n <- nrow(genotypes); p <- ncol(genotypes);

# Algo parms
seed <- 1; set.seed(seed)
kMax <- 5; tryNb <- 20; tol <- 1e-6
algoParms <- paste0('-seed', seed, '-try', tryNb, '-tol', tol)


################################################################################
# Fit mixtures
resName <- paste0(dataName, '-mixture-k', kMax)
resFile <- paste0(resName, algoParms, '.Rdata')

emList <- initList <- list(); par(mfrow=c(ceiling(sqrt(kMax)), round(sqrt(kMax))))
for(k in 1:kMax){
    cat('\n k =', k, ':')
    # Init tau
    initList[[k]] <- InitEM(genotypes=genotypes, k=k)
    if(k > 1){
      initList[[k]]$tau <- initList[[k]]$tau + 1/k
      initList[[k]]$tau <- initList[[k]]$tau / rowSums(initList[[k]]$tau)
    }
    emList[[k]] <- EMgenotypeMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                       k=k, tauInit=initList[[k]]$tau, tol=tol)
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
    apply(emList[[k]]$gamma, c(1, 2), sum, na.rm=TRUE)
    if(k > 1){for(tt in 1:tryNb){
      cat('\n k =', k, 'try', tt, ':')
      emTmp <- EMgenotypeMixtureHW(haplotypes=haplotypes, haploArray=haploArray,
                                   k=k, tauInit=rdirichlet(n, rep(1, k)), tol=tol)
      if(emTmp$logLik > emList[[k]]$logLik){cat('!!!'); emList[[k]] <- emTmp}
    }}
    plot(emList[[k]]$logLikPath, type='b', ylab='', xlab='', main=paste('k =', k))
    save(covar, pop, id, haplotypes, haploArray, genotypes, initList, emList, file=resFile)
}



