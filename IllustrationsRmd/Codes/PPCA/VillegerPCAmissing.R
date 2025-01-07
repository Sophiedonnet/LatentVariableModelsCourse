# PCA for Villeger data analysis: fish species in Terminos Lagoon (Gulf of Mexico)
# GrP13-Book-MultiAnalEcolData, chap 12

rm(list=ls()); par(pch=20); palette('R3')
seed <- 1; set.seed(seed)
library(mvtnorm);
source('../Functions/FunctionsProbabilisticPCA.R')
dirFig <- '../../Figures/'
dirData <- '../../Data/Villeger/'
exportFig <- TRUE
scale <- TRUE

# Data
dataName <- 'Villeger2012c'
traits <- read.table(paste0(dirData, dataName, '_AJ-traits.csv'), sep=';', dec=',', header=TRUE)
rownames(traits) <- traits[, 1]
Y <- as.matrix(traits[, -1]); n <- nrow(Y); p <- ncol(Y)
if(scale){Y <- scale(Y)}

# Missing values
missRate <- .2
R <- matrix(rbinom(n*p, 1, 1-missRate), n, p)
Yobs <- Y; Yobs[which(R==0)] <- NA

# pPCA
qMax <- 10; ppcaList <- list()
par(mfrow=c(ceiling(sqrt(qMax)), round(sqrt(qMax))))
for(q in 1:qMax){
  ppcaList[[q]] <- EMpPCAmiss(Yobs=Yobs, R=R, q=q)
  plot(ppcaList[[q]]$logLpath, main=q, ylim=quantile(ppcaList[[q]]$logLpath, probs=c(0.1, 1)), xlab='')
}

# Selection
par(mfrow=c(1, 1))
logL <- unlist(lapply(ppcaList, function(ppca){ppca$logL}))
aic <- unlist(lapply(ppcaList, function(ppca){ppca$aic}))
qAic <- which.max(aic)
bic <- unlist(lapply(ppcaList, function(ppca){ppca$bic}))
qBic <- which.max(bic)
plot(logL, type='b', pch=20, ylim=range(c(logL, aic, bic)))
lines(aic, type='b', pch=20, col=4); abline(v=qAic, col=4, lty=2)
lines(bic, type='b', pch=20, col=2); abline(v=qBic, col=2, lty=2)

# Imputation
q <- qBic; ppca <- ppcaList[[q]]; 
imputeAll <- ImputePPCAall(Yobs=Yobs, R=R, ppca=ppca)
Yfull <- Yobs; Yfull[which(R==0)] <- imputeAll$Yimp[which(R==0)]
plot(Y[which(R==0)], imputeAll$Yimp[which(R==0)], pch=20); abline(h=0, v=0, a=0, b=1)
for(h in 2:8){
  ppcaFull <- NaivePPCA(Y=Yfull, q=q)
  imputeEach <- ImputePPCAeach(Yfull=Yfull, R=R, ppca=ppcaFull)
  points(Y[which(R==0)], imputeEach$Yimp[which(R==0)], pch=20, col=h)
  Yfull <- Yobs; Yfull[which(R==0)] <- imputeEach$Yimp[which(R==0)]
}
