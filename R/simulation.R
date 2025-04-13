# simulation.R

simulate_negbin_data <- function(n, p, beta_true, r) {
  set.seed(123)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  eta <- X %*% beta_true
  mu <- exp(eta)
  y <- rnbinom(n, size = r, mu = mu)
  data.frame(y = y, X = X)
}
