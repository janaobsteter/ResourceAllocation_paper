library(mvtnorm)
library(ggplot2)

### Bivarite intensity
bivariateIntensity <- function(nTested, nSelected, rho) {
  #print(nTested)
  nNewBorn = 4320

  # Set bivariate distribution
  # Set parameters
  
  p1 = nTested / nNewBorn
  p2 = nSelected / nTested
  
  x = qnorm(1 - p1)
  y = qnorm(1 - p2*p1)
  #rho = 0.5 / 0.8
  
  lower = c(x, y)
  upper = c(Inf, Inf)
  corr = matrix(c(1, rho, rho, 1), nrow=2)
  cov = (0.5 / 0.8) * sqrt(1 * 1)
  sigma = matrix(c(1, cov, cov, 1), nrow=2)
  
  # Cumulative bivariate
  p = pmvnorm(mean = c(0,0), lower = lower, upper = upper, corr = corr) [1]
  p/(p1*p2)
}
### Bivarite intensity
bivariateIntensityBlank <- function(p1, p2, pT, cor) {
  #print(nTested)
  #nNewBorn = 4320

  # Set bivariate distribution
  # Set parameters
  
  
  x = qnorm(1 - p1)
  y = qnorm(1 - p2*p1)
  rho = cor
  
  lower = c(x, y)
  upper = c(Inf, Inf)
  corr = matrix(c(1, rho, rho, 1), nrow=2)
  cov = (0.5 / 0.8) * sqrt(1 * 1)
  sigma = matrix(c(1, cov, cov, 1), nrow=2)
  
  # Cumulative bivariate
  p = pmvnorm(mean = c(0,0), lower = lower, upper = upper, corr = corr) [1]
  p/pT
}

### Trivarite intensity
trivariateIntensity <- function(nTested, nMladi = 8, nSelected=5) {
  #print(nTested)
  nNewBorn = 4320

  # Set bivariate distribution
  # Set parameters
  
  p1 = nTested / nNewBorn
  p2 = nMladi / nTested
  p3 = nSelected / nMladi
  
  x = qnorm(1 - p1)
  y = qnorm(1 - p1*p2)
  z = qnorm(1 - p1*p2*p3)
  rho = 0.5 / 0.8
  
  lower = c(x, y)
  upper = c(Inf, Inf)
  corr = matrix(c(1, rho, rho, 1), nrow=2)
  cov = (0.5 / 0.8) * sqrt(1 * 1)
  sigma = matrix(c(1, cov, cov, 1), nrow=2)
  
  # Cumulative bivariate
  p = pmvnorm(mean = c(0,0), lower = lower, upper = upper, corr = corr) [1]
  p/(p1*p2)
}

computeZscores <- function(nTested) {
  nNewBorn = 4320
  nSelected = 5
  
  # Set bivariate distribution
  # Set parameters
  
  p1 = nTested / nNewBorn
  p2 = nSelected / nTested
  
  x = qnorm(1 - p1)
  y = qnorm(1 - p2)
  return(c(x,y))
}

par <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/ResourceAllocation/ParameterFile_Simulation.csv")
colnames(par) [4] <- "NoControl"
par$NoControl <- as.factor(par$NoControl)


par$IntensityTotal <- sapply(X = par$telMTotal, FUN = bivariateIntensity)
par$IntensityTotal <- round(par$IntensityTotal, 2)
par
ggplot(par, aes(x=NoControl, y=IntensityTotal, colour=G_P, group=G_P)) + geom_point(size = 2) + 
  geom_line(alpha=0.2)  #geom_bar(stat="identity", position="dodge")

zscore <- sapply(X = par$telMTotal, FUN=computeZscores)
par$zX <- zscore[1,]
par$zY <- zscore[2,]

ggplot(par, aes(x=NoControl, y=zX, colour=G_P, group=G_P)) + geom_line() #geom_bar(stat="identity", position="dodge")
ggplot(par, aes(x=NoControl, y=zY, colour=G_P, group=G_P)) + geom_line() #geom_bar(stat="identity", position="dodge")


# How changing the sd changed the correlation
V = matrix(data = c(1.0, 0.6, 0.6, 1.0), nrow = 2, byrow = TRUE)
cov2cor(V)
V = matrix(data = c(1.0, 0.6, 0.6, 4.0), nrow = 2, byrow = TRUE)
cov2cor(V)


cov = (0.5 / 0.9) * sqrt(1 * 2)
