# Probabilistic and phylogenetic PCA for the lizard data set

rm(list=ls()); palette('R3'); par(mfrow=c(1, 1), pch=20)
library(ade4); library(ape); library(phytools)
library(mvtnorm); library(fields)
source('../Functions/FunctionsProbabilisticPCA.R')

# Data
data("lizards"); traits <- lizards$traits; tree <- read.tree(text=lizards$hprA)
head(traits); n <- nrow(traits); p <- ncol(traits)
speciesName <- rownames(traits); traitsName <- colnames(traits)
plot.phylo(tree)
traits <- scale(traits); head(traits);
image.plot(1:p, 1:p, cov(traits), xlab='traits', ylab='traits', axes=0)
axis(1, at=1:p, labels=traitsName)
axis(2, at=1:p, labels=traitsName)
traitsDist <- as.matrix(dist(traits, diag=TRUE))

# PCA
pca <- prcomp(traits)
barplot(pca$sdev^2)

# PPCA
qMax <- p
ppcaList <- lapply(1:qMax, function(q){NaivePPCA(Y=traits, q=q, scale=FALSE)})
logL <- unlist(lapply(ppcaList, function(ppca){ppca$logL}))
bic <- unlist(lapply(ppcaList, function(ppca){ppca$bic}))
plot(logL, type='b', ylim=range(c(logL, bic)))
lines(bic, col=2, type='b')
qBic <- which.max(bic); abline(v=qBic, col=2)

qHat <- 2; # qHat <- qBic
ppca <- ppcaList[[qHat]]
# image.plot(1:p, 1:p, ppcaList[[qHat]]$Sigma, xlab='traits', ylab='traits', axes=0)
# axis(1, at=1:p, labels=traitsName)
# axis(2, at=1:p, labels=traitsName)
# plot(cov(traits), ppca$Sigma); abline(0, 1, h=0, v=0)
# plot(cor(traits), cov2cor(ppca$Sigma)); abline(0, 1, h=0, v=0)

# Cophenetic covariance
phyloDist <- cophenetic.phylo(tree)
plot(phyloDist, traitsDist)
phyloTime <- (max(phyloDist) - phyloDist)/2
phyloC <- phyloTime
phyloCinv <- solve(phyloC)
# cophCov <- vcvPhylo(tree, anc.nodes=FALSE)
# plot(cophCov, phyloTime); abline(0, 1)

# PhyloPCA
mu <- as.vector(1/(rep(1, n)%*%phyloCinv%*%rep(1, n))[1, 1] * (rep(1, n)%*%phyloCinv%*%traits))
phyloR <- 1/n * (t(traits-(rep(1, n)%o%mu)) %*% phyloCinv %*% (traits-(rep(1, n)%o%mu)))
plot(cor(traits), cov2cor(ppca$Sigma)); abline(0, 1, h=0, v=0)

# Removing phylo correlation
phyloCeig <- eigen(phyloC)
phyloCchol <- phyloCeig$vectors%*%diag(1/sqrt(phyloCeig$values))%*%t(phyloCeig$vectors)
traitsTilde <- phyloCchol%*%traits
traitsDistTilde <- as.matrix(dist(traitsTilde, diag=TRUE))
plot(as.data.frame(cbind(as.vector(phyloDist), as.vector(traitsDist), as.vector(traitsDistTilde))))

phyloPpcaList <- lapply(1:qMax, function(q){NaivePPCA(Y=traitsTilde, q=q, scale=FALSE)})
phyloLogL <- unlist(lapply(phyloPpcaList, function(ppca){ppca$logL}))
phyloLBic <- unlist(lapply(phyloPpcaList, function(ppca){ppca$bic}))
plot(phyloLogL, type='b', ylim=range(c(phyloLogL, phyloLBic)))
lines(phyloLBic, col=2, type='b')
qPhyloBic <- which.max(phyloLBic); abline(v=qPhyloBic, col=2)

