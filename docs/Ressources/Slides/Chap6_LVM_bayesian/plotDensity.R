library(ggplot2)
n = 1000
Z <- sample(c(0,1),n,replace=TRUE)
theta <- Z * rgamma(n,8,2) + (1-Z) * rgamma(n,10,1)
hist(theta)
df <- as.data.frame(theta)
names(df)<-'theta'



dens <- function(x){0.5 *dgamma(x,8,2) + 0.5*dgamma(x,10,1)}
y <- seq(0,20,length = 1000)
fy <- dens(y)
df2 <- as.data.frame(y)
names(df2)<-'y'
df2$fy <- fy


p <- ggplot(df, aes(x=theta)) +geom_histogram(aes(y=..density..), binwidth = 0.5,color="black", fill="red",alpha=0.1)
p <- p + geom_line(data=df2,aes(x = y,y=fy),color = 'red')
p
