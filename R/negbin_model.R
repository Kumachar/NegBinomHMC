#' Log-posterior for Negative Binomial Regression with Log-Overdispersion Parameter
#'
#' This function calculates the log-posterior for a Negative Binomial regression model,
#' where the overdispersion parameter is parameterized as phi = log(r).
#'
#' @param theta A numeric vector containing regression coefficients (beta) and the
#'   log-dispersion parameter (phi), i.e., theta = c(beta, phi).
#' @param X A numeric matrix of predictors.
#' @param y A numeric vector of response counts.
#' @param beta_mu A numeric vector of prior means for the regression coefficients.
#' @param beta_sigma A numeric vector or scalar representing the prior standard
#'   deviations for the regression coefficients.
#' @param alpha_r A numeric scalar representing the shape parameter of the gamma
#'   prior for r = exp(phi).
#' @param beta_r A numeric scalar representing the rate parameter of the gamma
#'   prior for r = exp(phi).
#'
#' @return A single numeric value representing the log-posterior.
#'
#' @export

log_posterior_negbin_log_r <- function(theta, X, y, beta_mu, beta_sigma, alpha_r, beta_r) {
  p <- ncol(X)
  beta <- theta[1:p]
  phi <- theta[p + 1]
  r <- exp(phi)
  
  eta <- X %*% beta
  mu <- exp(eta)
  
  # Log-likelihood
  log_likelihood <- sum(lgamma(y + r) - lgamma(r) - lgamma(y + 1) +
                          r * log(r) + y * eta - (r + y) * log(r + mu))
  
  # Log-prior for beta (normal)
  log_prior_beta <- sum(dnorm(beta, beta_mu, beta_sigma, log = TRUE))
  
  # Log-prior for r (gamma) with Jacobian adjustment
  log_prior_r <- dgamma(r, shape = alpha_r, rate = beta_r, log = TRUE) + phi
  
  return(log_likelihood + log_prior_beta + log_prior_r)
}

#' Gradient of the Log-posterior for Negative Binomial Regression with Log-Overdispersion
#'
#' This function computes the gradient of the log-posterior with respect to the regression
#' coefficients (beta) and the log-dispersion parameter (phi = log(r)).
#'
#' @param theta A numeric vector containing regression coefficients (beta) and the
#'   log-dispersion parameter (phi), i.e., theta = c(beta, phi).
#' @param X A numeric matrix of predictors.
#' @param y A numeric vector of observed counts.
#' @param beta_mu A numeric vector of prior means for the regression coefficients.
#' @param beta_sigma A numeric vector or scalar representing the prior standard
#'   deviations for the regression coefficients.
#' @param alpha_r A numeric scalar representing the shape parameter of the gamma
#'   prior for r = exp(phi).
#' @param beta_r A numeric scalar representing the rate parameter of the gamma
#'   prior for r = exp(phi).
#'
#' @return A numeric vector representing the gradient of the log-posterior with
#'   respect to theta = c(beta, phi).
#'
#' @export
grad_log_posterior_negbin_log_r <- function(theta, X, y, beta_mu, beta_sigma, alpha_r, beta_r) {
  p <- ncol(X)
  beta <- theta[1:p]
  phi <- theta[p + 1]
  r <- exp(phi)
  
  eta <- X %*% beta
  mu <- exp(eta)
  
  # Gradient with respect to beta
  grad_ll_beta <- t(X) %*% (y - ((r + y) * mu) / (r + mu))
  grad_prior_beta <- -(beta - beta_mu) / beta_sigma^2
  
  # Gradient with respect to phi
  # Check for small r to prevent numerical instability
  if (r < 1e-6) {
    # Return a large negative gradient to push phi away from very negative values
    grad_phi <- -1e6
  } else {
    grad_ll_r <- sum(digamma(y + r) - digamma(r) + log(r) + 1 -
                       log(r + mu) - (r + y) / (r + mu))
    grad_prior_r <- (alpha_r - 1) / r - beta_r
    grad_ll_phi <- grad_ll_r * r
    grad_prior_phi <- ((alpha_r - 1) - beta_r * r) * r + 1
    grad_phi <- grad_ll_phi + grad_prior_phi
  }
  
  # Combine gradients
  grad <- c(grad_ll_beta + grad_prior_beta, grad_phi)
  
  return(as.vector(grad))
}
