rm(list=ls())
setwd("~/Dropbox/Multiplex/Ecologie/Exposes/2017-06-Inecol-Xalapa/Code_R_plots/")

A=read.table("bluthgen_et_al_2004.txt")

et=read.csv("bluthgen_2004.csv")
head(A)
et2=read.csv("bluth.csv",skip=5,header=FALSE)

rowSums(A)
colSums(A)

et[1,-1]
abondance_col=as.numeric(as.vector(unlist(et[4,5:55])))
abondance_ligne=as.numeric(as.vector(et2[,3]))
et[46,]
head(et)



#binariser
A=as.matrix(A)
A[A>0]=1
library(blockmodels)

modLBM=BM_bernoulli("LBM",A)
modLBM$estimate()
modLBM
#classif
k=which.max(modLBM$ICL)
rcl=apply(modLBM$memberships[[k]]$Z1,1,which.max)
ccl=apply(modLBM$memberships[[k]]$Z2,1,which.max)

#affichage
G=graph_from_incidence_matrix(A)
plot(G,layout=layout_as_bipartite,vertex.label=NA,vertex.size=10,vertex.color=c(rcl,ccl+max(rcl)))
# pdf("bluthgen.pdf")
# dev.off()


#version resumee
pi=modLBM$model_parameters[[k]]$pi
piS=pi
piS[piS<.1]=0
G=graph_from_incidence_matrix(pi,weighted=TRUE)
alpha=c(table(rcl),table(ccl))
plot(G,layout=layout_as_bipartite,vertex.label=NA,vertex.size=alpha/2,vertex.color=1:(max(rcl)+max(ccl)),edge.width=E(G))
# pdf("bluthgen_sum.pdf")
# dev.off()



#lbm exemple
n=7
m=5
A=matrix(rbinom(n*m,1,.3),n,m)
class=c("orange","green","red","orange","red","red","green","yellow","black","yellow","black","black")
A[2,5]=1
G=graph_from_incidence_matrix(A)
plot(G,layout=layout.bipartite,vertex.color=class,vertex.label=NA)
# pdf("LBM_exemple.pdf")
# dev.off()
