where_save_file <- '/home/sophie/Dropbox/WORK_DROPBOX/ENSEIGNEMENT/2022-Orsay-MathsSV/Slides_LVM_CoursComplet/Chap2_LVM_Mixtures/plotsChap2'

library(mvtnorm)
library(plotly)
library(MASS)
library(ggplot2)

prop = c(1/2,1/3,1/6) 
mu = matrix(0,2,3)
mu[,2]=c(2,4)
mu[,3]=c(-2,6)

n = 100000 
Z = sample(1:3,n,prob=prop,replace=TRUE)
Sigma = matrix(0,2,2)
diag(Sigma)=c(1,1)
Y <- t(mu[,Z]) +rmvnorm(n,c(0,0),Sigma) 


### 
mymix <- as.data.frame(Y[1:1000,])
names(mymix) = c('Y1','Y2')
mymix$Z <- Z[1:1000]


filename=paste0(where_save_file,'/bivariategaussian.png')
ggplot(mymix,aes(x=Y1,y=Y2))+geom_point(data=subset(mymix,Z==1))
ggsave(filename = filename)

filename=paste0(where_save_file,'/bivariategaussianmixture.png')
ggplot(mymix,aes(x=Y1,y=Y2))+geom_point()
ggsave(filename = filename)

filename=paste0(where_save_file,'/bivariategaussianmixtureLabelled.png')
ggplot(mymix,aes(x=Y1,y=Y2))+geom_point(aes(colour = factor(Z))) + theme(legend.position = "none")
ggsave(filename = filename)





# the new part:
den3d <- kde2d(Y[,1], Y[,2])
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()
