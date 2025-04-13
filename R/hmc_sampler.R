# hmc_sampler.R

hmc_sampler <- function(log_post, grad_log_post, initial_beta, epsilon, L, n_iter, 
                        X, y, r, beta_mu, beta_sigma) {
  beta_current <- initial_beta
  samples <- matrix(NA, nrow = n_iter, ncol = length(beta_current))
  acceptances <- 0
  
  for (i in 1:n_iter) {
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
    
    # Negate momentum for symmetry
    momentum_proposed <- -momentum_proposed
    
    # Calculate acceptance probability
    current_U <- -log_post(beta_current, X, y, r, beta_mu, beta_sigma)
    current_K <- sum(momentum_current^2) / 2
    proposed_U <- -log_post(beta_proposed, X, y, r, beta_mu, beta_sigma)
    proposed_K <- sum(momentum_proposed^2) / 2
    
    log_acceptance_ratio <- current_U + current_K - proposed_U - proposed_K
    
    if (log(runif(1)) < log_acceptance_ratio) {
      beta_current <- beta_proposed
      acceptances <- acceptances + 1
    }
    
    samples[i, ] <- beta_current
  }
  
  acceptance_rate <- acceptances / n_iter
  list(samples = samples, acceptance_rate = acceptance_rate)
}
