# PCA for lizard data analysis
# GrP13-Book-MultiAnalEcolData, chap 12

rm(list=ls()); par(pch=20)
library(mvtnorm); library(mclust)
source('Functions-ProbabilisticPCA.R')
dirFig <- '../../Figures/'
exportFig <- TRUE
scale <- TRUE

# Data
dataName <- 'Lizards'
library(ade4); library(ape)
data("lizards"); traits <- lizards$traits; tree <- read.tree(text=lizards$hprA)
Y <- as.matrix(traits)
colnames(Y)[5] <- 'hatch.M'
colnames(Y)[7] <- 'matur.A'
colnames(Y) <- gsub(".", '', colnames(Y), fixed=TRUE)
head(Y)

# Description
if(exportFig){png(paste0(dirFig, dataName, '-scatter.png'))}
plot(as.data.frame(Y), pch=20)
if(exportFig){dev.off()}
cor(Y)

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

