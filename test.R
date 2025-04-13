# main.R

source("R/negbin_model.R")
source("R/hmc_sampler.R")
source("R/simulation.R")

# Simulation settings
n <- 200
p <- 3
beta_true <- c(0.5, -0.8, 1.2)
r <- 2

data_sim <- simulate_negbin_data(n, p, beta_true, r)

y <- data_sim$y
X <- as.matrix(data_sim[, -1])

# Priors
beta_mu <- rep(0, p)
beta_sigma <- 10

# HMC sampling parameters
epsilon <- 0.01  # step size
L <- 20          # number of leapfrog steps
n_iter <- 4000   # number of iterations

initial_beta <- rep(0, p)

result_hmc <- hmc_sampler(log_posterior_negbin, grad_log_posterior_negbin, 
                          initial_beta, epsilon, L, n_iter, 
                          X, y, r, beta_mu, beta_sigma)

# Acceptance rate
cat("Acceptance Rate:", result_hmc$acceptance_rate, "\n")

# Posterior summaries
posterior_samples <- result_hmc$samples
colnames(posterior_samples) <- paste0("beta_", 1:p)
summary(posterior_samples)

# Trace plot
par(mfrow = c(p, 1))
for (j in 1:p) {
  plot(posterior_samples[, j], type = "l", 
       main = paste("Trace plot for beta", j), ylab = "Value", xlab = "Iteration")
  abline(h = beta_true[j], col = "red", lwd = 2)
}

