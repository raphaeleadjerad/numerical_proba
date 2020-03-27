# EM for Gaussian mixture 

# Adjerad R. 

# packages ----------------------------------------------
library(ggplot2)


# Parameters --------------------------------------------
theta = list(pi = c(.2,.2,.6), mu = c(2,-5,2), sig = c(2,1,.1))
# Mixture of Gaussians, K = 3

# Function that calculates density f_theta(x) -----------
dnormmix <- function(x,theta){
  dens <- rep(0,length(x))
  K <- length(theta$pi)
  for (k in 1:K){
    dens <- dens + theta$pi[k]*dnorm(x,theta$mu[k],theta$sig[k])
  }                                         
  return(dens)
}
dnormmixObs <- function(obs,theta){
  n <- length(theta$pi)
  f <- sapply(1:n, function(i) theta$pi[i]*dnorm(obs,
                                                 mean = theta$mu[i],
                                                 sd = theta$sig[i]))
  dens <- rowSums(f)
  return(dens)
}


# Function that generates Gaussian mixture --------------
rnormmix <- function(n,theta){
  rnormix_real <- function(theta){
    K <- length(theta$pi)
    i <- sample(1:K, prob = theta$pi)
    X <- rnorm(1, theta$mu[i], theta$sig[i])
    return(X)
  }
  obs <- replicate(n,rnormix_real(theta))
  return(obs)
}
# Function to return log likelihood of Gaussian mixture -
lvraisnorm <- function(param,obs){    
  logvrais <- sum(log(dnormmix(obs,param)))
  return(logvrais)
}

# Visualizing likelihood for mixture of size K = 2 ------
theta <- list(pi=c(1,1)/2,mu=c(0,2.5),sig=c(1,1))
obs <- rnormmix(1000,theta)
mu1 <- seq(-1,4,by=.02)
mu2 <- seq(-2,4,by=.02)
L <- matrix(NA,nrow=length(mu1),ncol=length(mu2))
param <- theta
for (k in 1:length(mu1)){
  for (l in 1:length(mu2)){
    param$mu <- c(mu1[k],mu2[l]) 
    L[k,l] <- lvraisnorm(param,obs)
  } 
}
image(mu1, mu2, L, col = grey(0:255/255))
contour(mu1, mu2, L, add=T, nlevels=50, col='red')

# EM Algorithm ------------------------------------------

# Function to initialize parameters ---------------------
init_para <- function(obs, K = 3, size = 5){
  mu <- replicate(K,mean(sample(obs,size,replace=T)))
  sig <- replicate(K,sd(sample(obs,size,replace=T)))   
  theta.init <- list(pi=rep(1/K,K),mu=mu,sig=sig)
  return(theta.init)
  
}

# Function update.theta ---------------------------------
update.theta <- function(obs, theta){
  K <- length(theta$pi)
  n <- length(obs)
  # Calcul des alpha
  alpha <- matrix(nrow = n, ncol = K)
  for (k in 1:K){
    for(i in 1:n){
      alpha[i,k] <- theta$pi[k]*dnorm(obs[i], 
                                      mean = theta$mu[k], 
                                      sd = theta$sig[k])/dnormmix(obs[i],theta)
    }
  }
  # Calcul des pi
  pi.new <- sapply(1:K, function(k) mean(alpha[,k]))
  # Calcul des mu
  beta <- matrix(nrow = n, ncol = K)
  for (k in 1:K){
    for(i in 1:n){
      beta[i,k] <- obs[i]*alpha[i,k]/pi.new[k]
    }
  }
  mu.new <- sapply(1:K, function(k) mean(beta[,k]))
  # calcul des sigma
  gamma <- matrix(nrow = n, ncol = K)
  for (k in 1:K){
    for (i in 1:n){
      gamma[i,k] <- alpha[i,k]*((obs[i] - mu.new[k])^2)/pi.new[k]
    }
  }
  sig.new <- sqrt(sapply(1:K, function(k) mean(gamma[,k])))
  theta.new <- list(pi=pi.new, mu=mu.new, sig=sig.new)
  return(theta.new)
}

# Example of update.theta -----------------------------------
theta = list(pi=c(3,2,5)/10, mu=c(4,-5,0), sig=c(1.3,1,1.5))
obs <- rnormmix(500,theta)
update.theta(obs,theta)  # Note here we take the true theta to update theta 

# Algorithm EM function --------------------------------

algoEM <- function(obs, init_method = init_para, R=200, epsilon=0.001){   
  theta.old <- init_para(obs)
  theta.init <- init_para(obs)
  crit_arret <- FALSE
  log.vrais.old <- lvraisnorm(theta.init,obs)
  it <- 0
  while (!crit_arret && (it < R)){ 
    theta.new <- update.theta(obs, theta.old)
    log.vrais.new <- lvraisnorm(theta.new,obs)
    crit_arret <- (abs((log.vrais.new - log.vrais.old)/log.vrais.old) < epsilon)
    log.vrais.old <- log.vrais.new
    theta.old <- theta.new
    it <- it + 1
  }
  resultat <- list(emv = theta.new, nb.iter = it)
  return(resultat)
}

# Example of use of algoEM ----------------------------
res <- algoEM(obs)
res
theta
theta.em <- algoEM(obs)$emv
# The fit is far from perfect

# plot ------------------------------------------------
ymax <- max(c(dnorm(obs,mean(obs),sd(obs)),dnormmix(obs,theta)))
hist(obs,freq=FALSE,breaks = 15,ylim=c(0,ymax+0.03))
curve(dnorm(x,mean(obs),sd(obs)),add=TRUE,col='green')
curve(dnormmix(x,theta),add=TRUE,col='red')
curve(dnormmix(x,theta.em),add=TRUE,col='black')
# We fail to identify the three modes here

# Select model for best EM 
num_sim = 100
EMsim <- function(obs,K = 3,num_sim=100){
  model_good <- vector("list",num_sim)
  for (l in 1:num_sim)
    model_good[[l]] <- algoEM(obs)$emv
  logvrais <- sapply(model_good, function(x) lvraisnorm(x,obs))
  best.run <- which.max(logvrais)    
  theta <- model_good[[best.run]]   
  return(theta)
}
theta.em <- EMsim(obs)
ymax <- max(c(dnorm(obs,mean(obs),sd(obs)),dnormmix(obs,theta.em)))
hist(obs,freq=FALSE,breaks = 15,ylim=c(0,ymax+0.03))
curve(dnormmix(x,theta),add=TRUE,col='blue')
curve(dnormmix(x,theta.em),add=TRUE,col='red')
# Now we manage the identify the three modes better