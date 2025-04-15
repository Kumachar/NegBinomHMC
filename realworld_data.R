library(haven)
library(NegBinomHMC)
dat <- read_dta("https://stats.idre.ucla.edu/stat/stata/dae/nb_data.dta")
dat <- within(dat, {
  prog <- factor(prog, levels = 1:3, labels = c("General", "Academic", "Vocational"))
  id <- factor(id)
})

head(dat)
X <- model.matrix(~ math + prog, data = dat)
y <- dat$daysabs
r <- 2
p <- ncol(X)
# Define the log-posterior function for the negative binomial model
# Fit a negative binomial regression model
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

cat("Estimated beta (mean):", colMeans(result$samples), "\n")
cat("Acceptance rate:", result$acceptance_rate, "\n")
cat("Final epsilon:", result$final_epsilon, "\n")

