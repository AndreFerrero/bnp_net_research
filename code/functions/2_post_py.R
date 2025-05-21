post_py_sample <- function(N_predict,
                           alpha,
                           sigma,
                           data_py,
                           N_trunc = 10000,
                           base_sampler = function() rnorm(1)) {
  require(MCMCpack)
  
  # 1) Data summary
  n      <- length(data_py$obs_labels)
  K      <- data_py$n_clusters
  counts <- data_py$counts
  labels <- data_py$labels
  
  # 2) Mix weight
  R <- rbeta(1, n - sigma*K, alpha + sigma*K)
  
  # 3) Empirical part S
  wS <- as.numeric(rdirichlet(1, counts - sigma))
  S  <- list(weights = wS, atoms = labels)
  
  # 4) Simulate Q via py_sample_lab
  Q      <- py_sample_lab(N = N_trunc,
                          alpha = alpha + sigma*K,
                          sigma = sigma,
                          base_sampler = base_sampler)
  Q_obs  <- Q$obs_labels       # length N_trunc
  
  # 5) Prepare label→cluster map
  all_labs    <- c(S$atoms, Q_obs)
  unique_labs <- unique(all_labs)
  lab2cl      <- setNames(seq_along(unique_labs), unique_labs)
  
  # 6) Posterior predictive draws, sequentially from Q_obs
  pred_samples  <- numeric(N_predict)
  pred_clusters <- integer(N_predict)
  q_counter     <- 1
  
  for (i in seq_len(N_predict)) {
    if (runif(1) < R) {
      # from S
      lab <- sample(S$atoms, 1, prob = S$weights)
    } else {
      # from Q, in the exact order pre‐simulated
      if (q_counter > length(Q_obs)) {
        stop("Ran out of pre‐simulated Q samples (N_trunc too small)!")
      }
      lab <- Q_obs[q_counter]
      q_counter <- q_counter + 1
    }
    pred_samples[i]  <- lab
    pred_clusters[i] <- lab2cl[as.character(lab)]
  }
  
  list(
    R              = R,
    S              = S,
    Q              = Q,
    pred_samples   = pred_samples,
    pred_clusters  = pred_clusters,
    label_map      = lab2cl
  )
}

set.seed(123)
obs_data <- py_sample_lab(N = 30, alpha = 5, sigma = 0.3)
post_res <- post_py_sample(
  N_predict = 2000,
  alpha     = 20,
  sigma     = 0,
  data_py   = obs_data,
  N_trunc   = 1e4
)

# Plot predicted labels
hist(post_res$pred_samples,
     breaks = 30,
     main = "Posterior Predictive Samples",
     xlab = "Sampled Value", col = "skyblue", border = "white")
abline(v = obs_data$labels[1:obs_data$n_clusters],
       col = "red", lty = 2)
