#---Gaussian Process Mixture Model EM Algorithm---#
# Used in Chapter 3 of thesis

# kernel hyperparameter updates -- for squared exponential cov function + nugget
GP.K.par.updates <- function(Y, S, Gamma, mu, psi){ 
  G <- ncol(Gamma)
  n <- nrow(Y)
  p <- ncol(Y)
  
  sig2 <- psi$sig2
  l <- psi$l
  tau2 <- psi$tau2
  
  # updates via Nelder-Mead
  for (g in 1:G){
    w <- Gamma[, g]
    mu_g <- mu[g, ]
    
    # function to be optimized, theta is log parameters
    negQ <- function(theta){
      sig2_g <- exp(2*theta[1])
      l_g <- exp(theta[2])
      tau2_g <- exp(2*theta[3])
      
      U_g <- GPvecchia::createU(S, covparms = c(0, 1, sig2_g, l_g), nuggets = tau2_g, covmodel = "esqe")
      
      Q <- 0
      for (i in 1:n){
        Yc <- Y[i, ] - mu_g
        ll_Ug <- GPvecchia:::vecchia_likelihood_U(Yc, U_g)
        Q <- Q + w[i] * ll_Ug
      }
      return(-Q)
    }
    
    theta0 <- c(0.5*log(sig2[g]), log(l[g]), 0.5*log(tau2[g]))
    opt <- optim(theta0, negQ, method = "Nelder-Mead", 
                 control=list(  trace=0, parscale=rep(0.1,3), maxit=200, reltol=sqrt(.Machine$double.eps)))
    sig2[g] <- exp(2*opt$par[1])
    l[g] <- exp(opt$par[2])
    tau2[g] <- exp(2*opt$par[3])
  }
  return(psi=list(sig2 = sig2, l = l, tau2 = tau2 ))
}

# M-Step
GP.Mstep <- function(Y, S, Gamma, psi){
  G <- ncol(Gamma)
  n <- nrow(Y)
  p <- ncol(Y)
  
  # gamma sum over i
  cons <- apply( Gamma, 2, sum )
  # avoiding zero div
  cons[cons < 1e-12] <- 1e-12
  
  # pi estimate
  pi.dist <- cons / n
  
  # mean estimate -- mean 0 remains
  mu <- matrix(0, G, p)
  
  #hyper-par updates by gradient ascent 
  psi <- GP.K.par.updates(Y, S, Gamma, mu, psi)
  
  return( list( pi.dist=pi.dist, mu=mu, psi=psi))
  
}

# E-Step
GP.Estep <- function(par, Y, S){
  
  # parameters from M step
  pi.dist <- par$pi.dist
  mu <- par$mu
  psi <- par$psi #contains sigma^2, l and tau
  
  # number of states and obs
  G <- length(pi.dist)
  n <- nrow(Y)
  p <- ncol(Y)
  
  dens <- matrix(0, n, G)
  
  for( g in 1:G ){
    # parameters
    pi_g <- pi.dist[g]
    mu_g <- mu[g, ]
    sig2_g <- psi$sig2[g]
    l_g <- psi$l[g]
    tau2_g <- psi$tau2[g]
    
    # computing U_g once for each component
    U_g <- GPvecchia::createU(S, covparms = c(0, 1, sig2_g, l_g), nuggets = tau2_g, covmodel = "esqe")
    
    # approximated log-likelihood
    for (i in 1:n){
      Yc <- Y[i,] 
      ll_Ug <- GPvecchia:::vecchia_likelihood_U(Yc, U_g)
      dens[i,g] <- log(pi_g) + ll_Ug
    }
  }
  
  # for numerical stability
  max.vals <- apply( dens, 1, max )
  shift <- dens - max.vals
  Gamma.unnorm <- exp(shift)
  denom <- rowSums(Gamma.unnorm)
  
  # compute the gamma values 
  Gamma.mat <- Gamma.unnorm / denom
  
  # compute the log-likelihood 
  loglik <- sum( log(denom) + max.vals )
  
  return( list( Gamma.mat=Gamma.mat, loglik=loglik ) )
}

# hyper parameter initialisation
optimise.hyperpars <- function(X, Y, tol = 1e-4, maxstep = 1000){
  #sigma^2 = exp(2*theta[1])
  #l = exp(theta[2])
  #tau^2 = exp(2*theta[3])
  #D = as.matrix(dist(X)^2)
  
  # squared exponential covariance function and partial derivatives
  Kern <- function(D, gam1, gam2)
  {
    exp(2*gam1 - exp(-2*gam2)*D)
  }
  
  dKernsig <- function(D, gam1, gam2)
  {
    2*exp(2*gam1 - exp(-2*gam2)*D)
  }
  
  dKernell <- function(D, gam1, gam2)
  {
    2 * D * exp(2*gam1 - 2*gam2 - exp(-2*gam2)*D )
  }
  
  n <- nrow(Y)
  p <- ncol(Y)
  d <- ncol(X)
  nsub <- 0.5*p
  ny <- min(10, n)
  
  eta <- 0.00001 # learning rate
  grad <- numeric(3)
  eps <- Inf  
  steps <- 0
  
  # location points subset
  idx <- sort(sample(1:p, nsub, replace = FALSE))
  Xsub <- X[idx, , drop = FALSE]
  Dsub <- as.matrix(dist(Xsub))^2
  
  # curve subset
  curs <- sample(n, ny, replace = FALSE)
  Ysub <- Y[curs, idx, drop = FALSE]
  Syy <- crossprod(Ysub)
  
  sdv <- sd(as.numeric(Ysub))
  
  theta <- c( log(0.9*sdv), log( 0.1 * sqrt(max(Dsub))) , log(0.1*sdv))
  theta.old <- theta
  loglik.old <- -Inf
  
  # gradient ascent
  while( eps > tol & steps < maxstep + 1 ){
    K <- Kern(Dsub, theta[1], theta[2]) + diag(exp(2*theta[3]), nsub)
    cholK <- chol(K)
    Kinv <- chol2inv(cholK)
    
    loglik <- -(ny*sum(log(diag(cholK)))) - 0.5 * sum(Kinv * Syy)
    
    # common factors
    A <- Kinv %*% Syy %*% Kinv - ny * Kinv
    
    # derivatives
    # sigma
    dKsig <- dKernsig(Dsub, theta[1], theta[2])
    grad[1] <- 0.5 * sum(A * dKsig)
    
    # length
    dKell <- dKernell(Dsub, theta[1], theta[2])
    grad[2] <- 0.5 * sum(A * dKell)
    
    # nugget
    dKtau2 <- diag(2*exp(2*theta[3]), nsub)
    grad[3] <- 0.5 * sum(A * dKtau2)
    
    theta <- theta.old + eta * grad
    
    eps1 <- sum( abs(theta - theta.old) )
    if (!is.finite(loglik.old) || !is.finite(loglik)) {
      eps2 <- Inf
    } else {
      eps2 <- abs((loglik - loglik.old) / loglik.old) 
    }
    eps <- max(eps1, eps2)
    
    if( eps < 1e-4 ) converged <- TRUE else converged <- FALSE
    
    theta.old <- theta
    loglik.old <- loglik
    
    steps <- steps + 1 
  }
  
  return( list(theta = theta, par = c( exp(2*theta[1]), exp(theta[2]), exp(2*theta[3])), eps = eps))
}

# full initialisation
initialise.GP.EM <- function(G, Y, X){
  n <- nrow(Y)
  p <- ncol(Y)
  
  mu <- matrix(0, G, p)
  theta_hat <- matrix(NA, n, 3)
  
  for (i in 1:n){
    Yi <- Y[i, , drop = FALSE]
    opt <- optimise.hyperpars(X, Yi, tol = 1e-4, maxstep = 1000)
    theta_hat[i, ] <- log(opt$par)
  }
  
  km <- kmeans(theta_hat, centers = G, nstart = 50)
  z <- km$cluster
  
  Gamma.mat <- matrix(0, n, G)
  for (i in 1:n){
    Gamma.mat[i, z[i]] <- 1
  }
  
  mu <- matrix(0, G, p)
  
  sig2 <- numeric(G)
  l <- numeric(G)
  tau2 <- numeric(G)
  
  for (g in 1:G){
    sig2[g] <- exp(km$centers[g, 1])
    l[g]    <- exp(km$centers[g, 2])
    tau2[g] <- exp(km$centers[g, 3])
  }
  psi <- list(sig2 = sig2, l = l, tau2 = tau2 )
  return(list(Gamma.mat = Gamma.mat, mu=mu, psi=psi))
}

# ------ EM algorithm ----- #
EM.GP <- function( Y, X, G, m, maxstep=1000, tol=1e-6 ){
  
  # ordering + NN with maxmin ordering
  S <- GPvecchia::vecchia_specify(X, m, ordering = "maxmin", cond.yz = "y")
  
  # initialisation E-step
  init <- initialise.GP.EM(G, Y, X)
  Gamma <- init$Gamma.mat
  psi <- init$psi
  
  # M step with initial Gamma matrix and parameters
  mstep.out <- GP.Mstep(Y, S, Gamma, psi)
  
  # initialise old value of log-likelihood
  old.loglik <- -Inf
  eps <- Inf
  steps <- 0
  loglik.store <- numeric( maxstep )
  
  while( eps > tol & steps < maxstep + 1 )
  {
    # run E-step
    estep.out <- GP.Estep( mstep.out, Y, S )
    # extract info for M step
    Gamma <- estep.out$Gamma.mat
    # extract the value of the loglik
    new.loglik <- estep.out$loglik
    
    # update the parameter values
    mstep.out <- GP.Mstep(Y, S, Gamma, mstep.out$psi)
    
    # checking for convergence
    eps <- abs( ( new.loglik - old.loglik ) / new.loglik ) # will always be positive
    old.loglik <- new.loglik
    
    # increment step counter
    steps <- steps + 1
    
    loglik.store[ steps ] <- old.loglik
    cat("\n Value of the log-likelihood is: loglik = ", old.loglik )
  }
  
  if( eps < tol ) converged <- TRUE else converged <- FALSE
  
  if( !converged ) warning( "Algorithm did not converge in maximum allowed steps.")
  
  cat("\n Value of the log-likelihood is: loglik = ", old.loglik )
  
  # return order-sorted parameters
  o <- order(mstep.out$psi$l)
  mstep.out$mu <- mstep.out$mu[o, , drop = FALSE]
  mstep.out$psi$sig2 <- mstep.out$psi$sig2[o]
  mstep.out$psi$l <- mstep.out$psi$l[o]
  mstep.out$psi$tau2 <- mstep.out$psi$tau2[o]
  mstep.out$pi.dist <- mstep.out$pi.dist[o]
  mstep.out <- list(pi.dist=mstep.out$pi.dist, psi=mstep.out$psi)
  estep.out$Gamma.mat <- estep.out$Gamma.mat[, o, drop = FALSE]
  
  return( list( par=mstep.out, Gamma.mat=estep.out$Gamma.mat, 
                converged=converged, steps=steps,
                loglik.trace = loglik.store[1:steps], loglik = old.loglik ) )
}

# run EM 10 times and selecting highest log-likelihood solution
EM.GP.restarts <- function( Y, X, G, m, maxstep=1000, tol=1e-6, restarts=10)
{
  
  best <- list( loglik=-Inf )
  
  for( r in 1:restarts )
  {
    cat("\n Restart", r, "of", restarts, "\n")
    Z <- EM.GP( Y, X, G, m, maxstep, tol)
    if( Z$loglik > best$loglik )
    {
      best <- Z
      cat("\n*\t New maximum likelihood found: loglik = ", best$loglik )
    }
  }
  
  cat("\n")
  
  return( best )
}
