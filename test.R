# main.R

library(NegBinomHMC)
?hmc_sampler


# Test the adaptive sampler
n <- 200
p <- 3
beta_true <- c(1, 0, -1)
r <- 2
data <- simulate_negbin_data(n, p, beta_true, r, seed = 123)
X <- data$X
y <- data$y

result <- hmc_sampler(
  log_post = log_posterior_negbin,
  grad_log_post = grad_log_posterior_negbin,
  initial_beta = rep(0, p),
  initial_epsilon = 0.01,
  L = 10,
  n_iter = 1000,
  n_warmup = 500,
  X = X,
  y = y,
  r = r,
  beta_mu = 0,
  beta_sigma = 1,
  target_accept = 0.6
)

# Results
cat("True beta:", beta_true, "\n")
cat("Estimated beta (mean):", colMeans(result$samples), "\n")
cat("Acceptance rate:", result$acceptance_rate, "\n")
cat("Final epsilon:", result$final_epsilon, "\n")


# Posterior summaries
posterior_samples <- result$samples
colnames(posterior_samples) <- paste0("beta_", 1:p)
summary(posterior_samples)

# Trace plot
par(mfrow = c(p, 1))
for (j in 1:p) {
  plot(posterior_samples[, j], type = "l", 
       main = paste("Trace plot for beta", j), ylab = "Value", xlab = "Iteration")
  abline(h = beta_true[j], col = "red", lwd = 2)
}

