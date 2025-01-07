
##################### génération d'un SBM binaire #####################
rSBM <- function(n,Q,alpha,pi) {
  
  require (igraph)
  
  Z <- t(rmultinom(n, size = 1, prob = alpha))
  Znum <- Z %*% c(1:Q)
  
  Yvec <- rbinom(n^2,1,Z %*% pi %*% t(Z))
  X <- matrix(Yvec,n)
  X <- X * lower.tri(X) + t(X * lower.tri(X))
  diag(X) <- 0
  
  g <- graph.adjacency(X, mode='undirected', weighted=TRUE)
  return(list(X=X,g=g,cl=as.numeric(Znum),Z=Z))
}



rLBM <- function(n,p,alpha,beta,pi) {
  ### dim pi : Q lignes R colonnes
  require (igraph)
  Q = length(alpha)
  R = length(beta)
  
  Z <- t(rmultinom(n, size = 1, prob = alpha))  ## n lignes Q colonnes
  Y <- t(rmultinom(p, size = 1, prob = beta))  ## p lignes R colonnes
  Znum <- Z %*% c(1:Q)   #### de taille 
  Ynum <- Y %*% c(1:R)
  
  X = matrix(0,n,p)
  for (i in 1:n){for (j in 1:p){
    X[i,j] = rbinom(1,1,pi[Znum[i],Ynum[j]])
  }}

  
  g <- graph_from_incidence_matrix(X)
  return(list(X=X,g=g,cl_row=as.numeric(Znum),cl_col=as.numeric(Ynum),Z=Z,Y=Y))
}


