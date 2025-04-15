#' Random-Walk Metropolis Sampler for Negative Binomial Regression
#'
#' This function implements a random-walk Metropolis algorithm for sampling from the
#' posterior distribution of regression coefficients (\code{beta}) in a negative binomial
#' regression model. The algorithm proposes new values by adding a normally distributed noise
#' to the current parameter values and accepts/rejects proposals based on the Metropolis acceptance rule.
#'
#' @param log_post A function to compute the log-posterior of the model. The function
#'   should accept parameters in the following order: \code{log_post(beta, X, y, r, beta_mu, beta_sigma)}.
#' @param init A numeric vector providing the initial values for the parameter \code{beta}.
#' @param iterations An integer specifying the total number of iterations to run the sampler.
#' @param proposal_sd A numeric scalar (or vector of the same length as \code{init}) specifying
#'   the standard deviation(s) of the Gaussian proposal distribution.
#' @param X A numeric matrix representing the design matrix (including an intercept, if applicable).
#' @param y A numeric vector containing the observed count data.
#' @param r A numeric value that represents the dispersion parameter in the negative binomial model.
#' @param beta_mu A numeric value or vector indicating the mean(s) for the normal prior on \code{beta}.
#' @param beta_sigma A numeric value or vector indicating the standard deviation(s) for the normal prior on \code{beta}.
#'
#' @return A matrix of posterior samples for \code{beta}, where each row corresponds to one iteration
#'   and each column corresponds to a parameter.
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   # Generate synthetic data (assuming simulate_negbin_data is defined)
#'   n <- 200
#'   p <- 3
#'   beta_true <- c(1, 0, -1)
#'   r <- 2
#'   data <- simulate_negbin_data(n, p, beta_true, r, seed = 123)
#'   X <- data$X
#'   y <- data$y
#'
#'   # Set initial beta, number of iterations, and proposal standard deviation
#'   init <- rep(0, p)
#'   iterations <- 1000
#'   proposal_sd <- 0.1
#'
#'   # Run the Random-Walk Metropolis sampler
#'   chain <- rw_metropolis(
#'     log_post = log_posterior_negbin,
#'     init = init,
#'     iterations = iterations,
#'     proposal_sd = proposal_sd,
#'     X = X,
#'     y = y,
#'     r = r,
#'     beta_mu = 0,
#'     beta_sigma = 1
#'   )
#'
#'   # Check the posterior mean for beta
#'   colMeans(chain)
#' }
#'
#' @export
rw_metropolis <- function(log_post, init, iterations, proposal_sd, X, y, r, beta_mu, beta_sigma) {
  # init: initial parameter vector (for beta)
  chain <- matrix(NA, nrow = iterations, ncol = length(init))
  chain[1, ] <- init
  
  # Evaluate the log-posterior at the initial point
  current_log_post <- log_post(init, X, y, r, beta_mu, beta_sigma)
  
  # Iterative sampling loop
  for (i in 2:iterations) {
    # Propose new state by adding normally distributed noise
    proposal <- chain[i - 1, ] + rnorm(length(init), mean = 0, sd = proposal_sd)
    proposal_log_post <- log_post(proposal, X, y, r, beta_mu, beta_sigma)
    
    # Calculate acceptance probability
    accept_prob <- exp(proposal_log_post - current_log_post)
    
    # Accept or reject
    if (runif(1) < accept_prob) {
      chain[i, ] <- proposal
      current_log_post <- proposal_log_post
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
  }
  
  return(chain)
}
