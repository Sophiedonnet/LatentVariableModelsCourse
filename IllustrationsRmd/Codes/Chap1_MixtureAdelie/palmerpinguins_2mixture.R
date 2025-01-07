rm(list=ls())
library(ggplot2)
library(parallel)
library(viridis)
library(dplyr)
library(xtable)
#----------------------------------------- 
source('Functions_EM.R')
#---------------------------------------- 



#=======================================
#-------------- Data ---------------- 
#==========================================
library(palmerpenguins)

gg <- ggplot(penguins,aes(x=bill_length_mm,y = after_stat(density))) 
gg <- gg + geom_histogram(color = "white",fill="cyan4")
#gg <- gg + geom_density(color='cyan4')
gg + labs(title = "",
          x = "",
          y = "")



##########################################################
#------------  EM for all parameters
##########################################################

y <- penguins$bill_length_mm
y <- y[!is.na(y)]
#---------------- RUNs  EM from various starting points

KM <- kmeans(y, 2)

params0 = list(list(mu=c(40,50),var=c(5,5),p=0.5),
               list(mu=c(20,50),var=c(5,5),p=0.5),
               list(mu=c(35,70),var=c(5,5),p=0.6),
               list(mu=c(50,40),var=c(10,10),p=0.4))
params0[[5]] = list(mu=c(40,50),var=c(1,1),p=0.5)
params0[[6]] = list(mu=KM$centers,var=c(3,3),p=0.5)
res_all_run_all <-do.call(rbind,lapply(1:length(params0),function(run){
  res_EM_run  <- EM_2Mixture(params0[[run]],y,tol=10^{-7},estim_only_mu = FALSE)
  niter <- nrow(res_EM_run)
  res_EM_run <- cbind(res_EM_run,rep(run,niter))
  res_EM_run <- cbind(res_EM_run, c(1:nrow(res_EM_run)))
  names(res_EM_run)[7]='Init'
  names(res_EM_run)[8]='numIter'
  
  return(res_EM_run)}))
#------------------------------------------------
tabl <- matrix(unlist(params0),ncol=5,byrow=TRUE)
tabl <- cbind(tabl,rep(0,length(params0)))
row.names(tabl) = 1:length(params0)
for (init in 1:length(params0)){
  w <- which(res_all_run_all$Init==init)
  tabl[init,6]<- max(res_all_run_all$loglik[w])
}
xtable(tabl)


res_all_run_all$Init <- as.factor(res_all_run_all$Init)
res_all_run_all <- res_all_run_all %>% filter(numIter>3)
gg <- ggplot(res_all_run_all,aes(x=numIter,y=loglik,colour = Init)) + geom_line(linewidth=1.3) + geom_point(size=1.5) 
gg<- gg + theme(axis.text=element_text(size=15),axis.title=element_text(size=15),legend.text = element_text(size=15))
gg  + theme(legend.title = element_text(size=15)) + 
  labs(y = "Log likelihood",
       x = "Iterations")



#scale_fill_hue()$palette(5)  
##########################################################
#------------  EM on only MU 
##########################################################
theta_true = list(mu = c(35, 45))
theta_true$var = c(9,9)
theta_true$p = 0.36
Zsim = sample(c(1,2),length(y),prob=c(theta_true$p,1-theta_true$p),replace=TRUE)
ysim = theta_true$mu[Zsim] +sqrt(theta_true$var[Zsim])*rnorm(length(y))
hist(ysim)
N = 600

res_like = data.frame(mu1=double(),
                      mu2 = double(), 
                      ll = double())
grid_mu1 = seq(20,60,len=N)
grid_mu2 = seq(20,60,len=N)
grid_params <- as.data.frame(matrix(0,N*N,2))
names(grid_params) <- c('mu1', 'mu2') 
grid_params[,1]<-  rep(grid_mu1,N)
grid_params[,2] <- rep(grid_mu2,each=N)
res_loglik <- mclapply(1:N^2,function(i){
  params.i = theta_true
  params.i$mu <- c(grid_params[i,1],grid_params[i,2])
  logLikelihood(params.i,y=ysim)},mc.cores = 4)
res_like = grid_params
res_like[,3]  = unlist(res_loglik)
names(res_like)[3] <- 'loglik'

wmax = which.max(res_like[,3])
grid_params[wmax,]
 
#res_all_ru
#---------------- RUNs  EM from various starting points on only MU 

theta0 = theta_true
params0 <- lapply(1:5,function(i){theta0})
params0[[1]]$mu = c(25,25)
params0[[2]]$mu = c(30,25)
params0[[3]]$mu = c(55,40)
params0[[4]]$mu = c(25,55)
params0[[5]]$mu = c(25,45)




res_all_run <-do.call(rbind,lapply(1:length(params0),function(run){
  res_EM_run  <- EM_2Mixture(params0[[run]],ysim,tol=10^{-5},estim_only_mu = TRUE)
  niter <- nrow(res_EM_run)
  res_EM_run <- cbind(res_EM_run,rep(run,niter))
  res_EM_run <- cbind(res_EM_run, 1:niter)
  names(res_EM_run)[7]='Init'
  names(res_EM_run)[8]='numIter'
  
  return(res_EM_run)}))



#------------ PLOT ------------------

final_points <- res_all_run %>%
  group_by(Init) %>%
  filter(numIter == max(numIter)) %>%
  filter(Init %in% c(1, 2, 4))

# Ajouter les log-vraisemblances finales au dataframe des points finaux
final_points <- final_points %>%
  mutate(loglik_final = format(round(loglik, 2), nsmall = 2)) # Formater les valeurs


res_all_run$Init <- as.factor(res_all_run$Init)
res_all_run$numIter <- as.factor(res_all_run$numIter)
true_mu <- as.data.frame(matrix(theta_true$mu,nrow=1))
names(true_mu)  =c('mu1','mu2')
gg <- ggplot(res_like,aes(x=mu1,y=mu2,z=loglik)) + geom_tile(aes(fill=loglik)) + scale_fill_viridis() + stat_contour(color="white", size=0.25)
gg <- gg + geom_point(data=res_all_run, aes(x=mu1, y=mu2,z=NULL,colour=Init),size=3) + geom_line(data=res_all_run, aes(x=mu1, y=mu2,z=NULL,col=Init),linewidth=1.3) + scale_shape_manual(values=c(15, 20))
gg <- gg + geom_point(data=true_mu, aes(x=mu1, y=mu2,z=NULL),size=4,color='black') 
#gg <- gg+  geom_line(linewidth=1.3) + geom_point(size=1.5) 
gg<- gg + theme(axis.text=element_text(size=15),axis.title=element_text(size=15),legend.text = element_text(size=10))
gg <- gg  + theme(legend.title = element_text(size=15))+ guides(shape="none") + labs(x=expression(mu[1]),y=expression(mu[2]))

# Filtrer les dernières itérations pour les runs 1, 2 et 4

# Ajouter les annotations au graphique
gg <- gg + geom_text(data = final_points, 
            aes(x = mu1, y = mu2, label = loglik_final), 
            color = "black", size = 5, vjust = -1)
gg



init_labels <-c(
  "1" = "Init 1 (25, 25)",
  "2" = "Init 2 (30, 25)",
  "3" = "Init 3 (55, 40)",
  "4" = "Init 4 (25, 55)",
  "5" = "Init 5 (25, 45)"
)

# Modifier la légende de Init en ajoutant les valeurs initiales
gg <- gg + 
  scale_colour_discrete(
    labels = init_labels, # Libellés personnalisés
    name = "Initialisation" # Titre de la légende
  )
gg







