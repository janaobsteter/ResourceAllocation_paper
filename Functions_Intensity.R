library(truncnorm)
library(mvtnorm)


# This is to obtain the intensity of selection for a given proportion (height of the ordinate)
univariateIntensity <- function (p) {
  (exp(-0.5 * qnorm(1 - p)**2) * 1/(sqrt(2 * pi)))/p
}

univariateIntensity2 <- function (p) {
  dnorm(qnorm(p, lower.tail=FALSE)) / p
}

zScoreForProportion <- function (p) {
  # We use qnorm to get the Z score (cdf) for a given proportion
  qnorm(1 - p)
}



### Bivarite intensity
bivariateIntensity <- function(nTested, nSelected) {
  nNewBorn = 4320

  # Set bivariate distribution
  # Set parameters
  
  p1 = nTested / nNewBorn
  p2 = nSelected / nTested
  
  x = qnorm(1 - p1)
  y = qnorm(1 - p2)
  rho = 0.5 / 0.8
  
  lower = c(x, y)
  upper = c(Inf, Inf)
  corr = matrix(c(1, rho, rho, 1), nrow=2)
  
  # Cumulative bivariate
  p = pmvnorm(lower = lower, upper = upper, corr = corr) [1]
  p / (p1*p2)
}

bivariateIntensity(8,5)

# This is the equation of 
bivariateDistribution <- function (x, y, rho) {
  1/(2 * pi * sqrt(1 - rho**2)) * exp(-(1 / 2 * ( 1 - rho**2)) * (x**2 + y**2 + 2*rho*x*y)) / (p1 * p2)
}



### Plotting the distrivutions
library(GA)

#input
mu1<-1 #mean of X_1
mu2<-1 #mean of X_2
sigma11<-1 #variance of X_1
sigma22<- 1 #variance of X_2
sigma12< rho #covariance of X_1 and X_2 

#plot
x1 <- seq(mu1-3, mu1+3, length= 500)
x2 <- seq(mu2-3, mu2+3, length= 500)
z <- function(x1,x2){ z <- exp(-(sigma22*(x1-mu1)^2+sigma11*(x2-mu2)^2-2*sigma12*(x1-mu1)*(x2-mu2))/(2*(sigma11*sigma22-sigma12^2)))/(2*pi*sqrt(sigma11*sigma22-sigma12^2)) }
f <- outer(x1,x2,z)
persp3D(x1, x2, f, theta = 30, phi = 30, expand = 0.5)



library(plotly)

#input
mu1<-0 #mean of X_1
mu2<-0 #mean of X_2
sigma11 <- 1#variance of X_1
sigma22 <- 1 #variance of X_2
sigma12 < rho #covariance of X_1 and X_2 

#plot
x1 <- seq(mu1-3, mu1+3, length= 100)
x2 <- seq(mu2-3, mu2+3, length= 100)
z <- function(x1,x2){ z <- exp(-(sigma22*(x1-mu1)^2+sigma11*(x2-mu2)^2-2*sigma12*(x1-mu1)*(x2-mu2))/(2*(sigma11*sigma22-sigma12^2)))/(2*pi*sqrt(sigma11*sigma22-sigma12^2)) }
f <- t(outer(x1,x2,z))


plot_ly(x=x1,y=x2,z=f,type = "contour")%>%layout(xaxis=list(title="x1"),yaxis=list(title="x2"))



library(ggplot2)

#input
mu1<-2 #mean of X_1
mu2<-3 #mean of X_2
sigma11<-2 #variance of X_1
sigma22<-0.5 #variance of X_2
sigma12<-0.4 #covariance of X_1 and X_2 
x2<-4 #the fixed value of X_2

#plot
x1<-seq(round(mu1+sigma11*(x2-mu2)/sigma22,1)-3,round(mu1+sigma11*(x2-mu2)/sigma22,1)+3,length=100)
Curve<-dnorm(x1, mean =mu1+sigma11*(x2-mu2)/sigma22, sd = sigma11-(sigma12^2)/sigma11)
NormalDistData<-data.frame(x1,Curve)
ggplot()+geom_line(data=NormalDistData,aes(x=x1,y=Curve))+ scale_x_continuous(breaks = seq(mu1+sigma11*(x2-mu2)/sigma22-3,mu1+sigma11*(x2-mu2)/sigma22+3,by=1))+ylab(paste("f(x1|x2=",as.character(x2),")"))+ggtitle(paste("Pdf of X1 given X2=",as.character(x2)))
  