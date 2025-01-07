#--------------------------------- 
#setwd("~/Dropbox/Multiplex/Ecologie/Exposes/2017-06-Inecol-Xalapa/Code_R_plots")
#setwd("D:/Dropbox/Multiplex/Ecologie/Exposes/2017-06-Inecol-Xalapa/Code_R_plots")
#--------------------------------- 
library(ggplot2)
library(igraph)
source('func_sampling.R')


List_type_graph <- c("Bipartite","Nested")

epsilon <- 0.09
rho= .6
n <- 60
p <- 30; 
M_Nested = matrix(epsilon,4,3)
M_Nested[,1] = 0.8; 
M_Nested[1:3,2] =0.7 ; 
M_Nested[1:2,3] = 0.9; 



M = matrix(epsilon,4,3)
M[1,1] <- M[1,3] <- rho
M[2,2] <-  rho
M[3,1:2] <- M[4,3] <-  rho




for (type_graph in List_type_graph){
    print(type_graph)
    alpha <- switch(type_graph,"Bipartite"   = c(1/4,1/4,1/4,1/4), "Nested"      = c(1/4,1/4,1/4,1/4))
    beta <- switch(type_graph, "Bipartite"   = c(1/3,1/6,1/2),"Nested"      = c(1/3,1/6,1/2))

    pi <- switch(type_graph,"Bipartite"   = M, "Nested"      = M_Nested)


    lbm <- rLBM(n,p,alpha,beta,pi)
    save(lbm,alpha,beta,pi,file=paste('plots_lbm/data_',type_graph,'.Rdata'))
    
    R <- lbm$g
    layout <- layout_as_bipartite(R)
  
    par(mfrow=c(1,1))
    plot(R,layout=layout,vertex.size = 7,vertex.label=NA)
    file_save = paste('plots_lbm/',type_graph,'_graphe_without_colors.png',sep='')
    dev.copy(png,file_save)
    dev.off()
  
  
    
    
    plot(R,layout=layout,vertex.color=V(lbm$g)$type+1+c(lbm$cl_row,max(lbm$cl_row)+lbm$cl_col),vertex.size=7,vertex.label=NA)
    file_save = paste('plots_lbm/',type_graph,'_graphe_with_colors.png',sep='')
    dev.copy(png,file_save)
    dev.off()
  
  
  
  #graphe resume
    piS=pi
    piS[pi<0.1]=0
    G = graph_from_incidence_matrix(piS,weighted = TRUE)
    layout <- layout_as_bipartite(G)
    par(mfrow=c(1,1))
    plot(G,layout=layout,vertex.color=c(1:max(lbm$cl_row),max(lbm$cl_row)+1:max(lbm$cl_col)), vertex.size =c(alpha,beta)*100,vertex.label=NA)
    file_save = paste('plots_lbm/',type_graph,'_graphe_resume.png',sep='')
    dev.copy(png,file_save)
    dev.off()
  
  ######### ADJCENCY MATRIX

  Adj=lbm$X
  index_row = rep(1:dim(Adj)[1],each=dim(Adj)[2])
  index_col = rep(1:dim(Adj)[2],dim(Adj)[1])

  melted_Adj= data.frame(plants=index_row , animals= index_col)
  link = rep(-10,dim(Adj)[2]*dim(Adj)[1])
  for (k in 1:(dim(Adj)[2]*dim(Adj)[1])){link[k] = Adj[index_row[k],index_col[k]]}
  melted_Adj$link = link
  melted_Adj$index_row = n-index_row
  melted_Adj$index_col = index_col


  g<- ggplot(data = melted_Adj, aes(y=index_row, x=index_col, fill=link)) + geom_tile()
  g
  file_save = paste('plots_lbm/',type_graph,'_adja.png',sep='')
  ggsave(file_save, width = 20, height = 20, units = "cm")

  ########  REORDERED ADJACENCY MATRIX

  r = order(lbm$cl_row)
  Z_ord = sort(lbm$cl_row)
  
  s = order(lbm$cl_col)
  Y_ord = sort(lbm$cl_col)
  
  sep_row =n-which(diff(Z_ord)!=0)-0.5
  sep_col =which(diff(Y_ord)!=0)+0.5
  
  Adj_ord = lbm$X[r,s]

  melted_Adj_ord = data.frame(plants=index_row,animals =index_col)
  link_ord = rep(-10,dim(Adj)[2]*dim(Adj)[1])
  for (k in 1:(dim(Adj)[2]*dim(Adj)[1])){link_ord[k] = Adj_ord[index_row[k],index_col[k]]}
  melted_Adj_ord$link = link_ord
  melted_Adj_ord$index_row=n-index_row
  melted_Adj_ord$index_col =index_col



  g<- ggplot(data = melted_Adj_ord, aes(x=index_col, y=index_row, fill=link)) + geom_tile() +ggtitle("Reordered incidence matrix")
  g <- g + geom_hline(yintercept=sep_row,linetype="dashed",color='orange')+ geom_vline(xintercept=sep_col,linetype="dashed",color='orange')
  g<- g + theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
  g <- g + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
  g
  file_save = paste('plots_lbm/',type_graph,'_reordered_adja_with_groups.png',sep='')
  ggsave(file_save, width = 20, height = 20, units = "cm")




  g2<- ggplot(data = melted_Adj_ord, aes(x=index_col, y=index_row, fill=link)) + geom_tile() +ggtitle("Reordered indidence matrix")
  g2 <- g2 + theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
  g2 <- g2 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
  g2
  file_save = paste('plots_lbm/',type_graph,'_reordered_adja_without_groups.png',sep='')
  ggsave(file_save, width = 20, height = 20, units = "cm")

}
