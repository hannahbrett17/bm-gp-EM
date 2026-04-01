#---Brownian Motion Mixture Model EM Algorithm---#
# Used in Chapter 2 of thesis

# M-Step
BMMM.Mstep <- function( Gamma, Y, dt )
{
  G <- ncol(Gamma)
  n <- nrow(Y) # number of obs
  m <- ncol(Y) # number of increments of each BM
  
  # gamma sum over i
  cons <- apply( Gamma, 2, sum )
  # avoiding zero div
  cons[cons < 1e-12] <- 1e-12
  
  # pi estimate
  pi.dist <- cons / n
  
  # mean estimate
  mu <- numeric(G)
  for (l in 1:G){
    mu[l] <- sum(Gamma[,l] * Y) / (dt * m * cons[l])
  }
  
  # variance estimate
  sig2 <- numeric(G)
  for (l in 1:G) {
    diff <- (Y - dt * mu[l])^2
    sig2[l] <- sum(Gamma[,l] * diff) / (dt * m * cons[l])
  }
  
  
  return( list( pi.dist=pi.dist, psi=list( mean=mu, sig2=sig2 ) ))
  
}

# E-Step
BMMM.Estep <- function( par, Y, dt )
{
  
  # parameters computed in the M step
  pi.dist <- par$pi.dist
  psi <- par$psi
  
  # number of states and obs
  G <- length(pi.dist)
  n <- nrow(Y)
  m <- ncol(Y)
  
  # evaluate the density of the observations at all states, mean + sd scaled by dt 
  dens <- matrix(0, nrow = n, ncol = G)
  for (k in 1:G) {
    dens[,k] <- (log(par$pi.dist[k]) 
                 + rowSums(dnorm(Y, mean = dt*par$psi$mean[k], 
                                 sd = sqrt(dt*par$psi$sig2[k]), 
                                 log = TRUE))) 
  }
  
  
  # for numerical stability
  max.vals <- apply( dens, 1, max )
  shift <- dens - max.vals
  Gamma.unnorm <- exp(shift)
  denom <- rowSums(Gamma.unnorm)
  
  # --- compute the gamma values --- #
  Gamma.mat <- Gamma.unnorm / denom
  
  # --- compute the log-likelihood --- #
  loglik <- sum( log(denom) + max.vals )
  
  return( list( Gamma.mat=Gamma.mat, loglik=loglik ) )
}

# initialisation functions
# mean-based initialisation (section 2.4.2)
initialise.means <- function( G, Y, dt )
{
  m <- ncol(Y)
  n <- nrow(Y)
  
  # kmeans 
  km <- kmeans( rowSums(Y), centers = G, nstart = 100)
  # mean estimate
  mu <- km$centers / (dt * m) 
  
  z <- km$cluster
  sig2 <- numeric(G)
  pi.dist <- numeric(G)
  for( g in 1:G ) 
  {
    # full BM in cluster k
    Yg <- Y[z == g, , drop = FALSE]
    # variance
    sig2[g] <- var(as.numeric(Yg)) / dt 
    # pi
    pi.dist[g] <- nrow(Yg) / n
  }
  
  return( list( pi.dist=pi.dist, psi=list( mean=mu, sig2=sig2 ) ) )
}

# variance-based initialisation (section 2.4.3)
initialise.variance <- function( G, Y, dt )
{
  m <- ncol(Y)
  n <- nrow(Y)
  
  # kmeans
  km <- kmeans(log(rowMeans(Y^2)), centers = G, nstart = 100)
  
  z <- km$cluster
  mu <- numeric(G)
  sig2 <- numeric(G)
  pi.dist <- numeric(G)
  for( g in 1:G ) 
  {
    # full BM in cluster k
    Yg <- Y[z == g, , drop = FALSE]
    # mean
    mu[g] <- mean(Yg) / dt
    # variance
    sig2[g] <- var(as.numeric(Yg)) / dt 
    # pi
    pi.dist[g] <- nrow(Yg) / n
  }
  
  return( list( pi.dist=pi.dist, psi=list( mean=mu, sig2=sig2 ) ) )
}

# ------ EM algorithm ----- #
# Running EM algorithm
EM.BMMM <- function( Y, G, dt, maxstep=1000, tol=1e-6, init_func)
{
  
  Y <- Y[,,drop=FALSE]
  
  # initialise the model parameters
  mstep.out <- init_func( G, Y, dt )
  
  # initialise old value of log-likelihood
  old.loglik <- -Inf
  
  eps <- Inf
  
  steps <- 0
  
  loglik.store <- numeric( maxstep )
  
  while( eps > tol && steps < maxstep + 1 )
  {
    # run E-step
    estep.out <- BMMM.Estep( mstep.out, Y, dt )
    # extract info for M step
    Gamma <- estep.out$Gamma.mat
    # extract the value of the loglik
    new.loglik <- estep.out$loglik
    
    # update the parameter values
    mstep.out <- BMMM.Mstep( Gamma, Y, dt )
    
    # use a relative tolerance threshold
    eps <- abs( ( new.loglik - old.loglik ) / new.loglik + 1e-8) # will always be positive
    
    old.loglik <- new.loglik
    
    # increment step counter
    steps <- steps + 1
    
    loglik.store[ steps ] <- old.loglik
    
  }
  
  if( eps < tol ) converged <- TRUE else converged <- FALSE
  
  if( !converged ) warning( "Algorithm did not converge in maximum allowed steps.")
  
  cat("\n Value of the log-likelihood is: loglik = ", old.loglik )
  
  # return order-sorted means
  o <- order( mstep.out$psi$sig2 )
  mstep.out$psi$mean <- mstep.out$psi$mean[o]
  mstep.out$psi$sig2 <- mstep.out$psi$sig2[o]
  mstep.out$pi.dist <- mstep.out$pi.dist[o]
  estep.out$Gamma.mat <- estep.out$Gamma.mat[, o]
  
  return( list( par=mstep.out, Gamma.mat=estep.out$Gamma.mat, 
                converged=converged, steps=steps,
                loglik.trace = loglik.store[1:steps], loglik = old.loglik ) )
}


