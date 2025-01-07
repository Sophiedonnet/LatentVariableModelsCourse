# Illustration of the PLN model on Barents data

rm(list=ls()); par(lwd=2, pch=20, mfrow=c(1, 1))
library(PLNmodels); library(fields); library(RColorBrewer)
# See http://www.sthda.com/english/wiki/colors-in-r for RColorBrewer
exportFig <- TRUE

#-------------------------------------------------------------------------------
# Dir
dataName <- 'BarentsFish'
dataDir <- '../../Data/Barents/'
figDir <- '../../Figures/'
alpha <- .05
palette <- "PiYG"

# Data & fit
load(paste0(dataDir, dataName, '.Rdata'))
load(paste0(dataDir, dataName, '-All.Rdata'))
load(paste0(dataDir, dataName, '-modelList.Rdata'))
load(paste0(dataDir, dataName, '-plnList.Rdata'))
n <- nrow(Y); p <- ncol(Y); d <- ncol(X)

#-------------------------------------------------------------------------------
# Full model: covariates'effect
#-------------------------------------------------------------------------------
beta <- matrix(as.numeric(plnFull$model_par$B), d+1, p)
varBeta <- matrix(diag(vcov(plnFull)), d+1, p); 
rownames(beta) <- rownames(varBeta) <- colnames(modelList[[length(modelList)]]$X)
colnames(beta) <- colnames(varBeta) <- colnames(Y)
# image.plot(1:d, 1:p, beta[-1, ], xlab='covariates', ylab='species', cex.lab=1.5)

betaBrk <- seq(-8, 8, by=2)
# speciesOrder <- order(prcomp(t(beta[-1, ]))$x[, 1])
# heatmap(t(beta[-1, speciesOrder]), Rowv=NA, Colv=NA, 
#         breaks=betaBrk, col=brewer.pal(n=length(betaBrk)-1, name = palette))

if(exportFig){png(paste0(figDir, dataName, '-plnFull-beta.png'))}
image.plot(1:d, 1:p, beta[-1, ], xlab='', ylab='', axes=0, 
           breaks=betaBrk, col=brewer.pal(n=length(betaBrk)-1, name = palette))
axis(1, labels=colnames(Xscale), at=1:d, cex.axis=1.25)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

# Pseudo-stats
stat <- (beta / sqrt(varBeta))[-1, ] # Remove the intercept for the tests
maxBrk <- 6; statBrk <- seq(-maxBrk, maxBrk, by=2)
statPlot <- stat; 
statPlot[which(abs(stat) > maxBrk)] <- maxBrk*sign(stat[which(abs(stat) > maxBrk)])
statBrkLegends <- statBrk; statBrkLegends[1] <- '-Inf'; statBrkLegends[length(statBrk)] <- '+Inf'

if(exportFig){png(paste0(figDir, dataName, '-plnFull-stat.png'))}
image.plot(1:d, 1:p, statPlot, xlab='', ylab='', axes=0, 
           breaks=statBrk, lab.breaks=statBrkLegends, col=brewer.pal(n=length(statBrk)-1, name = palette))
axis(1, labels=colnames(Xscale), at=1:d, cex.axis=1.25)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

# Pseudo-p values
pVal <- 2*pnorm(abs(stat), lower.tail=FALSE)
pAdj <- matrix(p.adjust(pVal, method="BH"), d, p)
signif <- sign(stat) * (pAdj < alpha)

if(exportFig){png(paste0(figDir, dataName, '-plnFull-signif.png'))}
image.plot(1:d, 1:p, statPlot, xlab='', ylab='', axes=0, 
           breaks=statBrk, lab.breaks=statBrkLegends, col=brewer.pal(n=length(statBrk)-1, name = palette))
axis(1, labels=colnames(Xscale), at=1:d, cex.axis=1.25)
axis(2, labels=colnames(Y), at=1:p, las=2)
for(c in 1:d){for(j in 1:p){
  if(signif[c, j]==1){points(c+.1, j, pch=24, cex=1.5, bg=3)}
  if(signif[c, j]==-1){points(c-.1, j, pch=25, cex=1.5, bg=2)}
}}
if(exportFig){dev.off()}

#-------------------------------------------------------------------------------
# Model selection
#-------------------------------------------------------------------------------
dList <- unlist(lapply(modelList, function(model){ncol(model$X)}))
dimList <- p*(p+1)/2 + p*dList
penList <- round(dimList * log(n)/2)
elboList <- round(unlist(lapply(plnList, function(fit){fit$criteria$loglik})))
bicList <- round(unlist(lapply(plnList, function(fit){fit$criteria$BIC})))
iclList <- round(unlist(lapply(plnList, function(fit){fit$criteria$ICL})))
entropyList <- round(unlist(lapply(plnList, function(fit){fit$entropy})))
bicOrder <- order(bicList, decreasing=TRUE)
bicBest <- bicOrder[1]
plnBest <- plnList[[bicBest]]

if(exportFig){png(paste0(figDir, dataName, '-pln-bic.png'))}
par(new=FALSE)
plot(elboList, type='b', ylim=range(elboList, bicList+elboList[1]-bicList[1]), 
     pch=20, lwd=2, xlab='models', ylab='elbo', cex.lab=1.5)
for(c in 1:(d+1)){
  text(mean(which(dList==c)), min(elboList)+.6*diff(range(elboList)), paste0('d=', c), cex=1.25)
}
abline(v=.5+which(diff(dList)>0), col=1, lwd=1)
par(new=TRUE)
plot(bicList, type='b', xlab='', ylab='', col=2, pch=20, lwd=2, axes=0, 
     ylim=range(elboList)+(max(bicList)-max(elboList)))
axis(side=4, col=2, col.ticks=2, col.axis=2)
par(new=FALSE)
abline(v=which.max(bicList), lty=2, col=2, lwd=2)
if(exportFig){dev.off()}

# Three best models
for(m in bicOrder[1:3]){
  cat(m, '& & ', colnames(Xscale)[modelList[[m]]$model[-1]], 
      '&', elboList[m], '&', dimList[m], '&', penList[m], '&', bicList[m], '\\\\ \n')
}
# Null model
m <- 1
cat(m, '& (null) &', colnames(Xscale)[modelList[[m]]$model[-1]], 
    '&', elboList[m], '&', dimList[m], '&', penList[m], '&', bicList[m], '\\\\ \n')
# Full model
m <- 16
cat(m, '& (full) &', colnames(Xscale)[modelList[[m]]$model[-1]], 
    '&', elboList[m], '&', dimList[m], '&', penList[m], '&', bicList[m], '\\\\ \n')

#-------------------------------------------------------------------------------
# Covariance structure
#-------------------------------------------------------------------------------
summary(diag(plnNull$model_par$Sigma))
summary(diag(plnFull$model_par$Sigma))
summary(diag(plnBest$model_par$Sigma))

corBrk <- seq(-1, 1, length.out=11)

if(exportFig){png(paste0(figDir, dataName, '-plnNull-cor.png'))}
image.plot(1:p ,1:p, cov2cor(plnNull$model_par$Sigma), xlab='', ylab='', 
           breaks=corBrk, col=tim.colors(length(corBrk)-1), axes=0)
axis(1, labels=colnames(Y), at=1:p, las=2)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

if(exportFig){png(paste0(figDir, dataName, '-plnFull-cor.png'))}
image.plot(1:p ,1:p, cov2cor(plnFull$model_par$Sigma), xlab='', ylab='', 
           breaks=corBrk, col=tim.colors(length(corBrk)-1), axes=0)
axis(1, labels=colnames(Y), at=1:p, las=2)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

# # Prediction
# XbetaAll <- cbind(rep(1, nrow(Y)), Xscale)%*%plnFull$model_par$B
# image.plot(cov2cor(plnNull$model_par$Sigma))
# image.plot(cor(XbetaAll))
# plot(cor(XbetaAll), cov2cor(plnNull$model_par$Sigma))

if(exportFig){png(paste0(figDir, dataName, '-plnBest-cor.png'))}
image.plot(1:p ,1:p, cov2cor(plnBest$model_par$Sigma), xlab='', ylab='', 
           breaks=corBrk, col=tim.colors(length(corBrk)-1), axes=0)
axis(1, labels=colnames(Y), at=1:p, las=2)
axis(2, labels=colnames(Y), at=1:p, las=2)
if(exportFig){dev.off()}

plot(plnFull$model_par$Sigma, plnBest$model_par$Sigma, ylim=(range(plnNull$model_par$Sigma)))
points(plnFull$model_par$Sigma, plnNull$model_par$Sigma, col=2)
abline(h=0, v=0, a=0, b=1)

