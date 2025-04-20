#' Hamiltonian Monte Carlo (HMC) Sampler
#'
#' This function implements the Hamiltonian Monte Carlo algorithm to sample from a posterior
#' distribution.
#'
#' @param log_post A function that computes the log-posterior density.
#' @param grad_log_post A function that computes the gradient of the log-posterior density.
#' @param initial_theta A numeric vector of initial parameter values, e.g., c(beta, phi).
#' @param initial_epsilon A positive numeric value for the initial step size.
#' @param L A positive integer representing the number of leapfrog steps.
#' @param n_iter A positive integer specifying the number of iterations (after warm-up).
#' @param n_warmup A non-negative integer specifying the number of warm-up iterations.
#' @param X A predictor matrix.
#' @param y A response vector.
#' @param beta_mu A numeric vector specifying the prior mean for beta.
#' @param beta_sigma A numeric vector or scalar specifying the prior standard deviation for beta.
#' @param alpha_r A numeric scalar for the shape parameter of the gamma prior for r.
#' @param beta_r A numeric scalar for the rate parameter of the gamma prior for r.
#' @param target_accept A numeric value in (0,1) for the target acceptance probability.
#' @param adapt_gamma A positive numeric value controlling the adaptation rate.
#' @param adapt_kappa A numeric value in (0,1) for smoothing step size updates.
#' @param adapt_t0 A positive numeric value for stabilization during adaptation.
#'
#' @return A list containing:
#' \describe{
#'   \item{samples}{A matrix of posterior samples after warm-up.}
#'   \item{acceptance_rate}{The overall acceptance rate.}
#'   \item{final_epsilon}{The final adapted step size.}
#' }
#'
#' @export
hmc_sampler <- function(log_post, grad_log_post, initial_theta, initial_epsilon, L,
                        n_iter, n_warmup = 0, X, y, beta_mu, beta_sigma,
                        alpha_r, beta_r, target_accept = 0.65, adapt_gamma = 0.05,
                        adapt_kappa = 0.75, adapt_t0 = 10) {
  if (initial_epsilon <= 0) stop("initial_epsilon must be positive")
  if (L < 1 || L != floor(L)) stop("L must be a positive integer")
  if (n_iter < 1 || n_iter != floor(n_iter)) stop("n_iter must be a positive integer")
  if (n_warmup < 0 || n_warmup != floor(n_warmup)) stop("n_warmup must be non-negative")
  if (length(initial_theta) != (ncol(X) + 1)) stop("initial_theta dimension mismatch")
  if (target_accept <= 0 || target_accept >= 1) stop("target_accept must be in (0,1)")
  
  theta_current <- initial_theta
  epsilon <- initial_epsilon
  samples <- matrix(NA, nrow = n_iter, ncol = length(theta_current))
  acceptances <- 0
  
  log_epsilon <- log(epsilon)
  log_epsilon_bar <- log_epsilon
  H_bar <- 0
  mu <- log(10 * initial_epsilon)
  
  for (i in 1:(n_warmup + n_iter)) {
    theta_proposed <- theta_current
    momentum_current <- rnorm(length(theta_current), 0, 1)
    momentum_proposed <- momentum_current
    
    # Leapfrog integration
    grad <- grad_log_post(theta_proposed, X, y, beta_mu, beta_sigma, alpha_r, beta_r)
    if (any(is.nan(grad))) {
      momentum_proposed <- momentum_current  # Skip update if gradient is NaN
    } else {
      momentum_proposed <- momentum_proposed + (epsilon / 2) * grad
    }
    
    for (l in 1:L) {
      theta_proposed <- theta_proposed + epsilon * momentum_proposed
      if (l != L) {
        grad <- grad_log_post(theta_proposed, X, y, beta_mu, beta_sigma, alpha_r, beta_r)
        if (any(is.nan(grad))) {
          momentum_proposed <- momentum_current  # Skip update
        } else {
          momentum_proposed <- momentum_proposed + epsilon * grad
        }
      }
    }
    
    grad <- grad_log_post(theta_proposed, X, y, beta_mu, beta_sigma, alpha_r, beta_r)
    if (any(is.nan(grad))) {
      momentum_proposed <- momentum_current
    } else {
      momentum_proposed <- momentum_proposed + (epsilon / 2) * grad
    }
    
    # Calculate acceptance probability
    current_U <- -log_post(theta_current, X, y, beta_mu, beta_sigma, alpha_r, beta_r)
    current_K <- sum(momentum_current^2) / 2
    proposed_U <- -log_post(theta_proposed, X, y, beta_mu, beta_sigma, alpha_r, beta_r)
    if (!is.finite(proposed_U) || any(is.nan(momentum_proposed))) {
      log_acceptance_ratio <- -Inf
      alpha <- 0
    } else {
      proposed_K <- sum(momentum_proposed^2) / 2
      log_acceptance_ratio <- current_U + current_K - proposed_U - proposed_K
      alpha <- min(1, exp(log_acceptance_ratio))
    }
    
    # Accept/reject
    if (is.finite(log_acceptance_ratio) && log(runif(1)) < log_acceptance_ratio) {
      theta_current <- theta_proposed
      acceptances <- acceptances + 1
    }
    
    # Adapt epsilon during warm-up
    if (i <= n_warmup) {
      H_bar <- (1 - 1 / (i + adapt_t0)) * H_bar + (1 / (i + adapt_t0)) * (target_accept - alpha)
      log_epsilon <- mu - (sqrt(i) / adapt_gamma) * H_bar
      log_epsilon_bar <- adapt_kappa * log_epsilon + (1 - adapt_kappa) * log_epsilon_bar
      epsilon <- exp(log_epsilon)
    } else {
      epsilon <- exp(log_epsilon_bar)
    }
    
    # Store samples after warm-up
    if (i > n_warmup) {
      samples[i - n_warmup, ] <- theta_current
    }
  }
  
  acceptance_rate <- acceptances / (n_warmup + n_iter)
  if (acceptance_rate < 0.2 || acceptance_rate > 0.95) {
    warning("Final acceptance rate ", round(acceptance_rate, 2), " may indicate poor tuning")
  }
  
  list(samples = samples, acceptance_rate = acceptance_rate, final_epsilon = epsilon)
}
