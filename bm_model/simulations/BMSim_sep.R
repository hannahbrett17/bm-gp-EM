# Simulation I: Well-separated Brownian motion mixture
# Corresponds to Section 2.4.2 of thesis

# load EM functions
source("../bm_model/BrownianMotionEM.R")

# --------SIMULATING BMMM--------- #
# Simulating Brownian motion mixture model to apply EM alg

# uniform time grid for each BM
# interval length
t <- 1

# increment length
dt <- 0.001

# number of increments
m <- floor( t / dt )

# number of BMs
n <- 1000

# states
G <- 3

# parameters
pi.dist <- c(1/3, 1/3, 1/3)
mu <- c(-3, 0, 3)
sig2 <- c(0.25, 0.36, 0.49)

# assigning each curve to a component
set.seed(1234)
z <- sample(1:G, size=n, replace=TRUE, prob=pi.dist)

# Matrix of BM increments, sum of row is one BM
Y <- matrix(0, nrow=n, ncol=m)

# making Y 
for (i in 1:n) {
  g <- z[i] 
  Y[i, ] <- rnorm(m, mean = mu[g] * dt, sd = sqrt(sig2[g] * dt))
}

# Plotting simulated Brownian paths 
X <- t(apply(Y[1:250,], 1, cumsum))
X <- cbind(rep(0, 250), X)  
colours <- c("#49d9de", "#c71e8f", "orange")
path_colours <- colours[z]

# plotting
limits <- abs(max(mu) * t) + 2 * sqrt(max(sig2) * t)
matplot(seq(0, t, by=dt), t(X), type="l", xlab="t", ylab="Y(t)",
        main="Simulated Brownian Paths",
        xlim=c(0, t), ylim=limits * c(-1, 1), 
        lty=1, col=path_colours)
legend("topleft", legend=paste("Component", sort(unique(z))),
       col=colours, lty=1, bty="n")

# Running the EM algorithm 
em <- EM.BMMM(Y, G, dt, maxstep=1000, tol=1e-6, init_func = initialise.means)

# Results
em$par

# parameter estimation errors
mean_est_error <- Metrics::rmse(mu, em$par$psi$mean)
var_est_error <- Metrics::rmse(sig2, em$par$psi$sig2)
print(mean_est_error)
print(var_est_error)

# confusion matrix of component labels and ARI
z_est <- max.col(em$Gamma.mat)
table(true = z, est = z_est)
ari <- mclust::adjustedRandIndex(z, z_est)
print(ari)

# recovering number of components with BIC
G_max <- 10
BIC <- numeric(G_max)

# looping through different numbers of components and returning one with min BIC
for (g in 2:G_max){
  em <- EM.BMMM(Y, g, dt, maxstep=1000, tol=1e-6, init_func = initialise.means)
  #BIC
  loglik <- em$loglik
  p <- 3*g - 1 #3 from each component, pi sums to 1, last is determined
  BIC[g] <- p*log(n) -2*loglik
}

# BIC plot
best_G <- which.min(BIC[-1]) + 1
cat("Best number of states according to BIC:", best_G, "\n")
plot(2:G_max, BIC[-1], type='b',
     main='BIC vs Number of Components', xlab='G', ylab='BIC')
abline(v=best_G, col='red', lty=2)



