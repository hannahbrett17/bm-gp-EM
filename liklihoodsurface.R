# Two-component Gaussian mixture observed-data log-likelihood surface
# plot in section 1.2.3

set.seed(123)

# parameters
G <- 2
mix_prob <- c(0.5, 0.5)
mu_true <- c(-1.25, 1.25)
sigma_true <- c(1, 1)

# simulating data
n <- 500
z <- sample(1:G, size = n, replace = TRUE, prob = mix_prob)
y <- rnorm(n, mean = mu_true[z], sd = sigma_true[z])

# likelihood
loglik_gmm <- function(mu1, mu2, x, pi1 = 0.5, sigma = 1) {
  dens <- pi1 * dnorm(x, mean = mu1, sd = sigma) +
    (1 - pi1) * dnorm(x, mean = mu2, sd = sigma)
  sum(log(dens))
}

# mean values grid 
mu1_grid <- seq(-4, 4, length.out = 120)
mu2_grid <- seq(-4, 4, length.out = 120)

ll_mat <- outer(
  mu1_grid,
  mu2_grid,
  Vectorize(function(a, b) loglik_gmm(a, b, y))
)

# shifted loglik
ll_mat <- ll_mat - max(ll_mat)

# plot
filled.contour(
  x = mu1_grid,
  y = mu2_grid,
  z = ll_mat,
  color.palette = function(n) hcl.colors(n, "Mako"),
  nlevels = 50,
  xlab = expression(mu[1]),
  ylab = expression(mu[2]),
  main = "Observed-data log-likelihood surface",
  
  plot.axes = {
    axis(1)
    axis(2)
    
    contour(
      x = mu1_grid,
      y = mu2_grid,
      z = ll_mat,
      levels = pretty(range(ll_mat), 50),
      add = TRUE,
      drawlabels = FALSE,
      col = "black",
      lwd = 1
    )
    
    abline(a = 0, b = 1, lty = 2, lwd = 2)
    
    points(mu_true[1], mu_true[2], pch = 19, cex = 1.4)
    points(mu_true[2], mu_true[1], pch = 19, cex = 1.4)
  }
)
