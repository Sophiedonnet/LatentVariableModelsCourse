setwd("~/Dropbox/Multiplex/Ecologie/Exposes/2017-06-Inecol-Xalapa/Code_R_plots")

library(igraph)
library(xtable)

##### Symetric
n=4
X=matrix(rbinom(n^2,1,.2),nrow=n,ncol=n)
diag(X)=0
X[lower.tri(X)]=X[upper.tri(X)]
isSymmetric(X)


xtable(X,digits = 1)

G=graph_from_adjacency_matrix(X,mode="undirected")
par(mar = c(2, 2, 2, 2))
plot(G)

##### Directed

n=4
X=matrix(rbinom(n^2,1,.4),nrow=n,ncol=n)
diag(X)=0
isSymmetric(X)


xtable(X,digits = 0,label = NULL)
G=graph_from_adjacency_matrix(X,mode="directed")
par(mar = c(2, 2, 2, 2))
plot(G,edge.curved  = 0.5)

### LBM

n=4
m=7
X=matrix(rbinom(n*m,1,.4),nrow=n,ncol=m)

xtable(X)

G=graph_from_incidence_matrix(X)
names.nodes <- c(paste('R',1:n,sep=''),paste('C',1:m,sep=''))
col.nodes <- c(rep('salmon2',n),rep('darkolivegreen3',m))
plot(G,layout=layout_as_bipartite,vertex.label=names.nodes, vertex.size = 23, vertex.color=col.nodes )

pdf("graphe_bipartite.pdf")
plot(G,layout=layout_as_bipartite)
dev.off()


