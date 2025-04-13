#' Hamiltonian Monte Carlo (HMC) Sampler
#'
#' This function implements the Hamiltonian Monte Carlo algorithm to sample from a posterior distribution. It uses dual averaging during the warm-up period to adapt the step size (epsilon) based on a target acceptance probability.
#'
#' @param log_post A function that computes the log-posterior density given the current parameter values and the data.
#' @param grad_log_post A function that computes the gradient of the log-posterior density.
#' @param initial_beta A numeric vector of initial parameter values.
#' @param initial_epsilon A positive numeric value for the initial step size.
#' @param L A positive integer representing the number of leapfrog steps in each iteration.
#' @param n_iter A positive integer specifying the number of iterations (after warm-up) for sampling.
#' @param n_warmup A non-negative integer specifying the number of warm-up iterations during which the step size is adapted.
#' @param X A predictor matrix used in the log-posterior and its gradient calculations.
#' @param y A response vector used in the log-posterior and its gradient calculations.
#' @param r Additional data or hyperparameter required by the log-posterior function.
#' @param beta_mu A numeric vector specifying the prior mean for the beta parameters.
#' @param beta_sigma A numeric vector or scalar specifying the prior standard deviation for the beta parameters.
#' @param target_accept A numeric value in (0,1) indicating the desired target acceptance probability for adaptation.
#' @param adapt_gamma A positive numeric value controlling the adaptation rate in dual averaging.
#' @param adapt_kappa A numeric value in (0,1) for smoothing the step size updates.
#' @param adapt_t0 A positive numeric value acting as a stabilization constant during adaptation.
#'
#' @return A list containing:
#' \describe{
#'   \item{samples}{A matrix of posterior samples collected after the warm-up period.}
#'   \item{acceptance_rate}{The overall acceptance rate across all iterations (warm-up + sampling).}
#'   \item{final_epsilon}{The final adapted step size after warm-up.}
#' }
#'
#' @examples
#' \dontrun{
#' # Define the log-posterior and its gradient functions
#' log_post <- function(beta, X, y, r, beta_mu, beta_sigma) {
#'   # Compute log-posterior
#' }
#'
#' grad_log_post <- function(beta, X, y, r, beta_mu, beta_sigma) {
#'   # Compute gradient of log-posterior
#' }
#'
#' # Data setup (example)
#' X <- matrix(rnorm(100 * 3), ncol = 3)
#' y <- rnorm(100)
#' r <- NULL
#' beta_mu <- rep(0, 3)
#' beta_sigma <- rep(1, 3)
#'
#' initial_beta <- rep(0, 3)
#' initial_epsilon <- 0.1
#' L <- 10
#' n_iter <- 1000
#' n_warmup <- 500
#'
#' # Run the HMC sampler
#' result <- hmc_sampler(log_post, grad_log_post, initial_beta, initial_epsilon,
#'                       L, n_iter, n_warmup, X, y, r, beta_mu, beta_sigma)
#'
#' # Examine results
#' samples <- result$samples
#' acceptance_rate <- result$acceptance_rate
#' }
#'
#' @export

hmc_sampler <- function(log_post, grad_log_post, initial_beta, initial_epsilon, L, 
                        n_iter, n_warmup, X, y, r, beta_mu, beta_sigma, 
                        target_accept = 0.65, adapt_gamma = 0.05, adapt_kappa = 0.75, 
                        adapt_t0 = 10) {
  # Input validation
  if (initial_epsilon <= 0) stop("initial_epsilon must be positive")
  if (L < 1 || L != floor(L)) stop("L must be a positive integer")
  if (n_iter < 1 || n_iter != floor(n_iter)) stop("n_iter must be a positive integer")
  if (n_warmup < 0 || n_warmup != floor(n_warmup)) stop("n_warmup must be non-negative")
  if (length(initial_beta) != ncol(X)) stop("initial_beta dimension mismatch")
  if (target_accept <= 0 || target_accept >= 1) stop("target_accept must be in (0,1)")
  
  beta_current <- initial_beta
  epsilon <- initial_epsilon
  samples <- matrix(NA, nrow = n_iter, ncol = length(beta_current))
  acceptances <- 0
  
  # Dual averaging variables
  log_epsilon <- log(epsilon)
  log_epsilon_bar <- log_epsilon
  H_bar <- 0
  mu <- log(10 * initial_epsilon) # Heuristic: bias toward larger epsilon initially
  
  for (i in 1:(n_warmup + n_iter)) {
    beta_proposed <- beta_current
    momentum_current <- rnorm(length(beta_current), 0, 1)
    momentum_proposed <- momentum_current
    
    # Leapfrog integration
    momentum_proposed <- momentum_proposed + (epsilon / 2) * 
      grad_log_post(beta_proposed, X, y, r, beta_mu, beta_sigma)
    
    for (l in 1:L) {
      beta_proposed <- beta_proposed + epsilon * momentum_proposed
      if (l != L) {
        momentum_proposed <- momentum_proposed + epsilon * 
          grad_log_post(beta_proposed, X, y, r, beta_mu, beta_sigma)
      }
    }
    
    momentum_proposed <- momentum_proposed + (epsilon / 2) * 
      grad_log_post(beta_proposed, X, y, r, beta_mu, beta_sigma)
    
    # Calculate acceptance probability
    current_U <- -log_post(beta_current, X, y, r, beta_mu, beta_sigma)
    current_K <- sum(momentum_current^2) / 2
    proposed_U <- -log_post(beta_proposed, X, y, r, beta_mu, beta_sigma)
    if (!is.finite(proposed_U)) {
      log_acceptance_ratio <- -Inf
      alpha <- 0
    } else {
      proposed_K <- sum(momentum_proposed^2) / 2
      log_acceptance_ratio <- current_U + current_K - proposed_U - proposed_K
      alpha <- min(1, exp(log_acceptance_ratio))
    }
    
    # Accept/reject
    if (log(runif(1)) < log_acceptance_ratio) {
      beta_current <- beta_proposed
      acceptances <- acceptances + 1
    }
    
    # Adapt epsilon during warm-up
    if (i <= n_warmup) {
      H_bar <- (1 - 1 / (i + adapt_t0)) * H_bar + (1 / (i + adapt_t0)) * (target_accept - alpha)
      log_epsilon <- mu - (sqrt(i) / adapt_gamma) * H_bar
      log_epsilon_bar <- adapt_kappa * log_epsilon + (1 - adapt_kappa) * log_epsilon_bar
      epsilon <- exp(log_epsilon)
    } else {
      # After warm-up, use stabilized epsilon
      epsilon <- exp(log_epsilon_bar)
    }
    
    # Store samples after warm-up
    if (i > n_warmup) {
      samples[i - n_warmup, ] <- beta_current
    }
  }
  
  acceptance_rate <- acceptances / (n_warmup + n_iter)
  if (acceptance_rate < 0.2 || acceptance_rate > 0.95) {
    warning("Final acceptance rate ", round(acceptance_rate, 2), " may indicate poor tuning")
  }
  
  list(samples = samples, acceptance_rate = acceptance_rate, final_epsilon = epsilon)
}
