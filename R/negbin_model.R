# negbin_model.R

# Log-posterior for Negative Binomial regression
log_posterior_negbin <- function(beta, X, y, r, beta_mu, beta_sigma) {
  eta <- X %*% beta
  mu <- exp(eta)
  
  log_likelihood <- sum(lgamma(y + r) - lgamma(r) - lgamma(y + 1) +
                          r * log(r) + y * eta - (r + y) * log(r + mu))
  
  log_prior <- sum(dnorm(beta, beta_mu, beta_sigma, log = TRUE))
  
  return(log_likelihood + log_prior)
}

# Gradient of the log-posterior
grad_log_posterior_negbin <- function(beta, X, y, r, beta_mu, beta_sigma) {
  eta <- X %*% beta
  mu <- exp(eta)
  grad_ll <- t(X) %*% (y - ((r + y) * mu) / (r + mu))
  grad_lp <- -(beta - beta_mu) / beta_sigma^2
  
  return(as.vector(grad_ll + grad_lp))
}
