# PCA for climate data analysis
# GrP13-Book-MultiAnalEcolData, chap 12

rm(list=ls()); par(pch=20)
library(mvtnorm); library(mclust)
source('Functions-ProbabilisticPCA.R')

# Data
dataName <- 'Climate'
climate <- read.table('climate.csv', sep=';', header=TRUE, row.names=1)
Y <- as.matrix(climate)

# dataName <- 'Ohrazeni'
# load('../Ohrazeni/ohrazeni.rda')
# Y <- ohrazeni$traitdata[, -(1:3)]
# Y <- Y[-which(rowSums(is.na(Y))>0), ]

head(Y)
scale <- FALSE

# # Scaling
# climateMean <- colMeans(Y)
# climateSd <- apply(Y, 2, sd)
# if(scale){Y <- scale(Y)}

# # HC
# hc <- hclust(dist(scale(Y)), method='ward.D2')
# plot(hc)

# PCA
pca <- prcomp(Y, scale=scale)
plot(pca)
plot(pca$x, col=0)
text(pca$x, col=1, labels=rownames(Y))
abline(v=0, h=0)

# # GMM
# gmm <- Mclust(Y)
# gmm
# plot(gmm, what='BIC')
# plot(gmm, what='classification')

# pPCA
ppcaList <- list()
for(q in 1:ncol(Y)){
  ppcaList[[q]] <- NaivePPCA(Y, q=q, scale=scale)
}
logL <- unlist(lapply(ppcaList, function(ppca){ppca$logL}))
bic <- unlist(lapply(ppcaList, function(ppca){ppca$bic}))
aic <- unlist(lapply(ppcaList, function(ppca){ppca$aic}))
plot(logL, type='b', col=1, ylim=range(c(logL, bic, aic)))
points(bic, type='b', col=2)
points(aic, type='b', col=4)
abline(v=which.max(bic), col=2)
# barplot(eigen(ppcaList[[ncol(Y)]]$Sigma)$values)

q <- 5
fit1 <- NaivePPCA(Y, q=q, scale=FALSE)
fit2 <- EMpPCA(Y, q=q)
plot(fit1$Sigma, fit2$Sigma); abline(0, 1)
max(abs(fit1$Sigma - fit2$Sigma))
max(abs(fit1$B - fit2$B))
