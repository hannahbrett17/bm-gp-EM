# Restart sensitivity of GP EM
# Corresponds to section 3.6.1 of thesis

#load EM functions
source("GaussianProcessEM.R")

# --------SIMULATING GPMM--------- #
# Simulating GP mixture model to apply EM alg

# Covariance function
sqex_kernel <- function(X, sig2, l, tau2)
{
  D <- as.matrix(dist(X)^2)
  K <- sig2 * exp( -D / l^2 ) 
  diag(K) <- diag(K) + tau2
  return(K)
}

# Simulating GP mixture model

# number of curves
n <- 75

# number of time points
p <- 100

# inputs
X <- seq(-1,1,length.out=p)
X <- as.matrix(X)

# number of components
G <- 3

# parameters
pi.dist <- c(1/3, 1/3, 1/3)
sig2 <- c(1, 1.5, 2)
l <- c(0.1, 0.2, 0.4)
tau2 <- c(1e-4, 1e-4, 1e-4)
mu <- matrix(0, G, p)
true_par <- list(pi.dist=pi.dist, sig2=sig2, l=l, tau2=tau2)

# assigning each curve to a component
set.seed(123)
z <- sample(1:G, size=n, replace=TRUE, prob=pi.dist)

# making Y
Y <- matrix(0, n, p)
for (i in 1:n){
  g <- z[i]
  K_g <- sqex_kernel(X, sig2[g], l[g], tau2[g])
  Y[i, ] <- as.numeric(mvtnorm::rmvnorm(1, mean = mu[g, ], sigma = K_g))
}

#---Running the EM algorithm 10 times for m=5 neighbours and calculating metrics---#
em_fulldata <- EM.GP.restartsdata(Y, X, G, m, true_par, z, restarts =10)
em_fulldata$summary
