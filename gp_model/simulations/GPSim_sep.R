# Simulation I: Separated GP mixture
# Corresponds to Section 3.5.1 of thesis

# load EM functions
source("../GaussianProcessEM.R")

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
n <- 50

# number of time points
p <- 100

# inputs
X <- seq(-1,1,length.out=p)
X <- as.matrix(X)

# number of components
G <- 2

# parameters
pi.dist <- c(1/2, 1/2)
sig2 <- c(1, 2)
l <- c(0.1, 0.5)
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
  Y[i, ] <- as.numeric(mvtnorm::rmvnorm(1, mean = mu[g, ], sigma = K_g))
}

# plotting mixture model
colours <- c("dodgerblue", "red")
path_colours <- colours[z[31:40]]
matplot(X, t(Y[31:40,]), type = "l", lwd=2,lty = 1, xlab = "x", ylab = "Y(x)", col = path_colours, 
        main = "Simulated GPs")
legend("topright", legend=paste("Component", sort(unique(z[31:40]))),
       col=colours, lty=1, lwd=2, bty="n")

#---Running the EM algorithm for m=5 neighbours---#
m <- 5
em <- EM.GP.restarts( Y, X, G, m, maxstep=1000, tol=1e-6, restarts=10)
em$par

# parameter estimation errors
sig2_error <- Metrics::rmse(sig2, em$par$psi$sig2)
l_error <- Metrics::rmse(l, em$par$psi$l)
tau2_error <- Metrics::rmse(tau2, em$par$psi$tau2)
print(sig2_error)
print(l_error)
print(tau2_error)

# confusion matrix of component labels and ARI
z_est <- max.col(em$Gamma.mat)
table(true = z, est = z_est)
ari <- mclust::adjustedRandIndex(z, z_est)
print(ari)

