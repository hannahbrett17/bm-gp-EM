# Effects of different m on Vecchia approximation
# Corresponds to section 3.6.2 of thesis

# load EM functions
source("GaussianProcessEM.R")

# Covariance function
sqex_kernel <- function(X, sig2, l, tau2)
{
  D <- as.matrix(dist(X)^2)
  K <- sig2 * exp( -D / l^2 ) 
  diag(K) <- diag(K) + tau2
  return(K)
}

# input grid
# number of time points
p <- 100

# inputs
X <- seq(-1,1,length.out=p)
X <- as.matrix(X)


#---Comparing likelihood for different m---#
# parameters
sig2 <- 1
l <- 0.1
tau2 <- 1e-4
mu <- rep(0, p)

# true covariance matrix
K_true <- sqex_kernel(X, sig2, l, tau2)

# making Y, one observation only
set.seed(123)
Y_single <- as.numeric(mvtnorm::rmvnorm(1, mu, K_true))

# Exact log-likelihood
ll_exact <- mvtnorm::dmvnorm(Y_single, mu, K_true, log=TRUE)

# Vecchia approximation for increasing m
m_vals_ll <- c(1, 3, 5, 8, 10)
for (m in m_vals_ll){
  S <- GPvecchia::vecchia_specify(X, m, ordering="maxmin", cond.yz="y")
  
  U <- GPvecchia::createU(S, covparms=c(0, 1, sig2, l), nuggets=tau2, covmodel="esqe")
  
  ll_vec <- GPvecchia:::vecchia_likelihood_U(Y_single, U)
  
  cat(sprintf("m = %2d | Vecchia = %.4f | Exact = %.4f | diff = %.6f\n",
              m, ll_vec, ll_exact, abs(ll_vec - ll_exact)))
}


#---Comparing EM results for different m---#
# Simulating mixture model (same as sim II, section 3.5.2) for reduced size n=50

# number of curves
n <- 50

# number of components
G <- 2

# parameters
pi.dist <- c(1/2, 1/2)
sig2 <- c(1, 1.2)
l <- c(0.2, 0.25)
tau2 <- c(1e-4, 1e-4)
mu <- matrix(0, G, p)

# assigning each curve to a component
set.seed(123)
z <- sample(1:G, size=n, replace=TRUE, prob=pi.dist)

# making Y
Y <- matrix(0, n, p)
for (i in 1:n){
  g <- z[i]
  K_g <- sqex_kernel(X, sig2[g], l[g], tau2[g])
  Y[i, ] <- as.numeric(mvtnorm::rmvnorm(1, mu[g, ], K_g))
}

# Running EM for different values of m
m_vals_em <- c(1, 3, 5, 8)
results <- data.frame(m=integer(), loglik=numeric(), ari=numeric())
em <-  vector("list", length(m_vals_em))

for (j in 1:length(m_vals_em)){
  m <- m_vals_em[j]
  em_m <- EM.GP.restarts(Y, X, G, m=m, maxstep=500, tol=1e-5, restarts=5)
  
  z_est <- max.col(em_m$Gamma.mat)
  
  ari_m <- mclust::adjustedRandIndex(z, z_est)
  results <- rbind(results, data.frame(m=m, loglik=em_m$loglik, ari=ari_m))
  em[[j]] <- em_m
}
print(results)
