
library(ggplot2)
library(sna)
setwd("~/WORK/RECHERCHE/FORMATIONS_CHERCHEURS/FORMATIONS_RESEAUX/IntroductionToBlockModels/Slides_SBM_LBM/Code_R_plots_for_slides")

############" SIMULATION FOR SBM
source('func_sampling.R')

epsilon <- .2
n <- 50
Q <- 3
type_graph <- "Bipartite"

alpha <- switch(type_graph,
                "Affiliation" = c(1/5,3/5,1/5),
                "Affiliation2" = c(.2,.6,.2),
                "Bipartite"   = c(1/4,1/4,1/4,1/4),
                "Etoile"      = c(.15,.35,.15,.35))
pi <- switch(type_graph,
             "Affiliation" = matrix(c(1-epsilon,epsilon,epsilon,epsilon,1-epsilon,epsilon,epsilon,epsilon,1-epsilon),3,3),
             "Affiliation2" = matrix(c(epsilon,epsilon/10,epsilon/10,epsilon/10,epsilon,epsilon/10,epsilon/10,epsilon/10,epsilon),3,3),
             "Bipartite"   = matrix(c(epsilon,1-epsilon,epsilon,epsilon,1-epsilon,epsilon,epsilon,epsilon,epsilon,epsilon, epsilon,1-epsilon,epsilon,epsilon,1-epsilon,epsilon),4,4),
             "Etoile"      = matrix(c(1-epsilon,1-epsilon,epsilon,epsilon,1-epsilon,epsilon,epsilon,epsilon,epsilon,epsilon,1-epsilon,1-epsilon,epsilon,epsilon,1-epsilon,epsilon),4,4))


Q=length(alpha)

sbm <- rSBM(n,Q,alpha,pi)
plot(graph_from_adjacency_matrix(sbm$X,mode="undirected"))


library(reshape2)
melted_sbm <- melt(sbm$X)
head(melted_sbm)
ggplot(data = melted_sbm, aes(x=Var1, y=Var2, fill=value)) + geom_tile()


################# matrix reordered according to groups
w = order(sbm$cl)

sbm_ord = sbm$X[w,w];

dim(sbm_ord)
melted_sbm_ord <- melt(sbm_ord)
head(melted_sbm_ord)
ggplot(data = melted_sbm_ord, aes(x=Var1, y=Var2, fill=value)) + geom_tile()+theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                                                                                             size = 12, hjust = 1))

theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                 size = 12, hjust = 1))
