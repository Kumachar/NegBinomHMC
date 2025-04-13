#' Simulate Negative Binomial Data for Regression
#'
#' This function generates synthetic data for Negative Binomial regression. It creates a predictor matrix
#' (including an intercept) and simulates count responses from a Negative Binomial distribution, using a
#' linear predictor based on the provided true regression coefficients.
#'
#' @param n A positive integer specifying the number of observations to simulate.
#' @param p A positive integer specifying the number of predictors (including the intercept). The length of \code{beta_true}
#' should be equal to \code{p}.
#' @param beta_true A numeric vector of true regression coefficients of length \code{p}. The first element is the intercept.
#' @param r A positive numeric value representing the dispersion parameter for the Negative Binomial distribution.
#' @param seed An optional integer to set the seed for reproducibility. If \code{NULL}, the random number generator is not seeded.
#'
#' @return A list containing:
#' \describe{
#'   \item{y}{A numeric vector of simulated count responses.}
#'   \item{X}{A numeric matrix of predictors including an intercept column.}
#' }
#'
#' @examples
#' \dontrun{
#' # Set simulation parameters
#' n <- 100       # Number of observations
#' p <- 3         # Number of predictors (including intercept)
#' beta_true <- c(0.5, -1.2, 0.8)  # True regression coefficients
#' r <- 1.5       # Dispersion parameter
#' seed <- 123    # Seed for reproducibility
#'
#' # Simulate data
#' sim_data <- simulate_negbin_data(n, p, beta_true, r, seed)
#' y <- sim_data$y
#' X <- sim_data$X
#'
#' # Display a preview of the generated data
#' head(X)
#' head(y)
#' }
#'
#' @export
simulate_negbin_data <- function(n, p, beta_true, r, seed = NULL) {
  if (n < 1 || n != floor(n)) stop("n must be a positive integer")
  if (p < 1 || p != floor(p)) stop("p must be a positive integer")
  if (r <= 0) stop("r must be positive")
  if (length(beta_true) != p) stop("beta_true must have length p")
  
  if (!is.null(seed)) set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  eta <- X %*% beta_true
  eta <- pmin(pmax(eta, -100), 100) # Prevent overflow
  mu <- exp(eta)
  y <- rnbinom(n, size = r, mu = mu)
  list(y = y, X = X)
}
