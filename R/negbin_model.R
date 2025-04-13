#' Log-posterior for Negative Binomial Regression
#'
#' This function calculates the log-posterior for a Negative Binomial regression model.
#' It combines the log-likelihood of the Negative Binomial distribution with a normal prior on
#' the regression coefficients.
#'
#' @param beta A numeric vector of regression coefficients.
#' @param X A numeric matrix of predictors.
#' @param y A numeric vector of response counts.
#' @param r A numeric value (or vector) representing the dispersion parameter for the Negative Binomial distribution.
#' @param beta_mu A numeric vector of prior means for the regression coefficients.
#' @param beta_sigma A numeric vector or scalar representing the prior standard deviations for the regression coefficients.
#'
#' @return A single numeric value representing the log-posterior.
#'
#' @examples
#' \dontrun{
#' # Example data setup
#' X <- matrix(rnorm(100), ncol = 1)
#' beta <- 0.5
#' y <- rnbinom(100, size = 1, mu = exp(X * beta))
#' r <- 1
#' beta_mu <- 0
#' beta_sigma <- 1
#'
#' # Compute the log-posterior
#' log_posterior_negbin(beta, X, y, r, beta_mu, beta_sigma)
#' }
#'
#' @export
log_posterior_negbin <- function(beta, X, y, r, beta_mu, beta_sigma) {
  eta <- X %*% beta
  mu <- exp(eta)
  
  log_likelihood <- sum(lgamma(y + r) - lgamma(r) - lgamma(y + 1) +
                          r * log(r) + y * eta - (r + y) * log(r + mu))
  
  log_prior <- sum(dnorm(beta, beta_mu, beta_sigma, log = TRUE))
  
  return(log_likelihood + log_prior)
}

#' Gradient of the Log-posterior for Negative Binomial Regression
#'
#' This function computes the gradient of the log-posterior with respect to the regression coefficients
#' for a Negative Binomial regression model. The gradient is useful for optimization routines and sampling methods
#' such as Hamiltonian Monte Carlo.
#'
#' @param beta A numeric vector of regression coefficients.
#' @param X A numeric matrix of predictors.
#' @param y A numeric vector of observed counts.
#' @param r A numeric value (or vector) representing the dispersion parameter for the Negative Binomial distribution.
#' @param beta_mu A numeric vector of prior means for the regression coefficients.
#' @param beta_sigma A numeric vector or scalar representing the prior standard deviations for the regression coefficients.
#'
#' @return A numeric vector representing the gradient of the log-posterior with respect to beta.
#'
#' @examples
#' \dontrun{
#' # Example data setup
#' X <- matrix(rnorm(100), ncol = 1)
#' beta <- 0.5
#' y <- rnbinom(100, size = 1, mu = exp(X * beta))
#' r <- 1
#' beta_mu <- 0
#' beta_sigma <- 1
#'
#' # Compute the gradient of the log-posterior
#' grad_log_posterior_negbin(beta, X, y, r, beta_mu, beta_sigma)
#' }
#'
#' @export
grad_log_posterior_negbin <- function(beta, X, y, r, beta_mu, beta_sigma) {
  eta <- X %*% beta
  mu <- exp(eta)
  grad_ll <- t(X) %*% (y - ((r + y) * mu) / (r + mu))
  grad_lp <- -(beta - beta_mu) / beta_sigma^2
  
  return(as.vector(grad_ll + grad_lp))
}
