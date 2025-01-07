## -------------------------------------------------------------------------------
library(sbm)
library(ggplot2)
library(stringr)
library(GGally)
library(bipartite)

#---------------------------------------
where_fig = "/home/sophie/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2024/2024_06_GTdataAgro/plots/Vacher"

## ----import dataset-------------------------------------------------------------
data("fungusTreeNetwork")
str(fungusTreeNetwork,  max.level = 1)

############""BIPARTITE####

#------------Plot Matrix Tree Fungis--------------------------
M <- fungusTreeNetwork$fungus_tree
rownames(M) <- str_sub(fungusTreeNetwork$fungus_names, 1, 30)
colnames(M) <- str_sub(fungusTreeNetwork$tree_names,1,30)
g <- plotMyMatrix(t(M), dimLabels = list(row = 'Trees', col = 'Fungis'), plotOptions = list(rowNames = TRUE,colNames = TRUE))
g + theme(axis.text=element_text(size=5))

#-------------------- Plot bipartite
plotweb(M,text.rot = 90)






#################################################
#------------------------------- BINARY
##################################################

## ----tree_tree_binary network---------------------------------------------------
tree_tree_binary <- 1 * (fungusTreeNetwork$tree_tree != 0)
row.names(tree_tree_binary) <- colnames(tree_tree_binary) <- str_sub(fungusTreeNetwork$tree_names, 1, 30)
g <- plotMyMatrix(tree_tree_binary , dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g + theme(axis.text=element_text(size=9))


tree_tree_binary_net = network(tree_tree_binary , directed = FALSE)
network.vertex.names(tree_tree_binary_net) = rownames(tree_tree_binary)
ggnet2(tree_tree_binary_net,  label = TRUE, label.size = 3)
## ----simpleSBM print,eval = TRUE, echo = TRUE-----------------------------------

mySimpleSBM <- tree_tree_binary %>%
  estimateSimpleSBM("bernoulli", 
                    dimLabels ='Trees', 
                    estimOptions = list(verbosity = 2, plot=FALSE))


## ----simpleSBMfit---------------------------------------------------------------

pBin <- mySimpleSBM$storedModels %>%  ggplot() + aes(x = nbBlocks, y = ICL)  + geom_line() + geom_point(alpha = 0.5) +theme(axis.text=element_text(size=9))
pBin <- pBin + geom_segment(aes(x = mySimpleSBM$nbBlocks, y = min(mySimpleSBM$storedModels$ICL) , xend = mySimpleSBM$nbBlocks, yend =  mySimpleSBM$ICL))
pBin  <- pBin + scale_x_discrete(limits=c("1", "2","3","4","5","6","7","8","9"))
ggsave(paste0(where_fig,"/FungusTree_Binary_ICL.png"),
  plot = pBin,
  scale = 1,
  width = 20,
  units = c("cm"))
## ----simpleSBMfit plot1---------------------------------------------------------
g <- plot(mySimpleSBM, type = "expected", dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=9))
ggsave(paste0(where_fig,"/FungusTree_Binary_expected.png"),
       plot = g,
       scale = 1,
       width = 20,
       units = c("cm"))

g <- plot(mySimpleSBM, type = "data", dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=9))
ggsave(paste0(where_fig,"/FungusTree_Binary_dataordered.png"),
       plot = g,
       scale = 1,
       width = 20,
       units = c("cm"))




#################################################
#------------------------------- Poisson
##################################################

## ----tree_tree network plot data------------------------------------------------

tree_tree_weighted <- fungusTreeNetwork$tree_tree
row.names(tree_tree_weighted) <- colnames(tree_tree_weighted) <- str_sub(fungusTreeNetwork$tree_names, 1, 30)
g <- plotMyMatrix(tree_tree_weighted , dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=7))

ggsave(paste0(where_fig,"/FungusTree_Poisson_data.png"),
       plot = g ,
       scale = 1,
       width = 20,
       units = c("cm"))
tree_tree_weighted_net = network(tree_tree_weighted , directed = FALSE)
network.vertex.names(tree_tree_weighted_net) = rownames(tree_tree_binary)
ggnet2(tree_tree_weighted_net,  label = TRUE, label.size = 3)

## ----simpleSBM Poisson, eval   = TRUE, echo = TRUE------------------------------
mySimpleSBMPoisson <- tree_tree_weighted  %>%
  estimateSimpleSBM("poisson", directed = FALSE,
                    estimOptions = list(verbosity = 0 , plot = FALSE),
                    dimLabels = c('Trees'))


## ----simpleSBMfitPoisson plot1--------------------------------------------------
pPoisson <- mySimpleSBMPoisson$storedModels %>%  ggplot() + aes(x = nbBlocks, y = ICL)  + geom_line() + geom_point(alpha = 0.5) +theme(axis.text=element_text(size=9))
pPoisson <- pPoisson + geom_segment(aes(x = mySimpleSBMPoisson$nbBlocks, y = min(mySimpleSBMPoisson$storedModels$ICL) , xend = mySimpleSBMPoisson$nbBlocks, yend =  mySimpleSBMPoisson$ICL))
pPoisson <- pPoisson + scale_x_discrete(limits=c("1", "2","3","4","5","6","7","8","9"))
pPoisson
ggsave(paste0(where_fig,"/FungusTree_Poisson_ICL.png"),
       plot = pPoisson,
       scale = 1,
       width = 20,
       units = c("cm"))

## ----simpleSBMfit plot1---------------------------------------------------------
g <- plot(mySimpleSBMPoisson, type = "expected", dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g + theme(axis.text=element_text(size=7))
ggsave(paste0(where_fig,"/FungusTree_Poisson_expected.png"),
       plot = g,
       scale = 1,
       width = 20,
       units = c("cm"))

g <- plot(mySimpleSBMPoisson, type = "data", dimLabels = list(row = 'Trees', col = 'Trees'), plotOptions = list(legend = TRUE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=7))
ggsave(paste0(where_fig,"/FungusTree_Poisson_dataordered.png"),
       plot = g,
       scale = 1,
       width = 20,
      units = c("cm"))

#################################################
#------------------------------- Poisson with Covar
##################################################
mySimpleSBMCov<- estimateSimpleSBM(
  netMat = as.matrix(tree_tree),
  model = 'poisson',
  directed =FALSE,
  dimLabels =c('Trees'), 
  covariates  = fungusTreeNetwork$covar_tree,
  estimOptions = list(verbosity = 0))


## ----select SBM covar, echo=TRUE, eval = TRUE-----------------------------------
mySimpleSBMCov$nbBlocks


## ----extract param SBM poisson covar, echo=TRUE, eval = TRUE--------------------
mySimpleSBMCov$connnectParam
mySimpleSBMCov$blockProp
mySimpleSBMCov$memberships
mySimpleSBMCov$covarParam


#################################################
#------------------------------- Bipartite
##################################################
## ----plot incidence-------------------------------------------------------------

fungus_tree <- t(fungusTreeNetwork$fungus_tree)
row.names(fungus_tree)  <- str_sub(fungusTreeNetwork$tree_names, 1, 30)
colnames(fungus_tree)  <- str_sub(fungusTreeNetwork$fungus_names, 1, 30)
pBip <- plotMyMatrix(fungus_tree, dimLabels=list(col = 'Fungis',col = 'Trees'), plotOptions = list(legend = FALSE,rowNames = TRUE,colNames = TRUE))
pBip <- pBip + theme(axis.text.x = element_text( size = 5),
             axis.text.y = element_text( size = 5))  
ggsave(paste0(where_fig,"/FungusTree_Biparite_data.png"),
       plot = pBip ,
       scale = 1,
       width = 30,
       units = c("cm"))


## ----tree_fungi_bipartite network, eval=TRUE, echo = TRUE-----------------------
myBipartiteSBM <- estimateBipartiteSBM(
  netMat =fungus_tree,
  model = 'bernoulli',
  dimLabels=c(col = 'Fungis',row = 'Trees'),
  estimOptions = list(verbosity = 0, plot=FALSE))




## ----plot bipartite estim-------------------------------------------------------
plot(myBipartiteSBM, type = "data")
g <- plot(myBipartiteSBM, type = "expected", dimLabels = list(col = 'Fungis',row = 'Trees'), plotOptions = list(legend = FALSE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=5))
ggsave(paste0(where_fig,"/FungusTree_Bipartite_expected.png"),
       plot = g,
       scale = 1,
       width = 30,
       units = c("cm"))

## ----plot bipartite expect------------------------------------------------------
g <- plot(myBipartiteSBM, type = "data", dimLabels = list(col = 'Fungis',row = 'Trees'), plotOptions = list(legend = FALSE,rowNames = TRUE,colNames = TRUE))
g <- g + theme(axis.text=element_text(size=5))
ggsave(paste0(where_fig,"/FungusTree_Bipartite_dataordered.png"),
       plot = g,
       scale = 1,
       width = 30,
       units = c("cm"))

## ----plot BM network tree fungi,  echo=TRUE, eval = TRUE------------------------
plot(myBipartiteSBM, type = "meso", plotOptions = list(vertex.size=c(0.5,1) , edge.width  = 1))

