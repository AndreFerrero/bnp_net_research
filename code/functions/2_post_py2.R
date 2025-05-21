post_py_sample <- function(N_predict,
                           alpha,
                           sigma,
                           data_py,
                           N_trunc = 10000,
                           base_sampler = function() rnorm(1)) {
  require(MCMCpack)
  
  # Step 1: Extract empirical summary
  n        <- length(data_py$obs_labels)
  K        <- data_py$n_clusters
  counts   <- data_py$counts
  labels   <- data_py$labels[1:K]
  
  # Step 2: Sample mixing weight
  R <- rbeta(1,
             shape1 = n - sigma * K,
             shape2 = alpha + sigma * K)
  
  # Step 3: Sample empirical component (S)
  dir_weights <- as.numeric(rdirichlet(1, counts - sigma))
  S <- list(
    weights = dir_weights,
    atoms   = labels
  )
  
  # Step 4: Sample new Pitman-Yor tail (Q)
  Q <- py_sample_lab(
    N            = N_trunc,
    alpha        = alpha + sigma * K,
    sigma        = sigma,
    base_sampler = base_sampler
  )
  
  # Step 5: Create a label-to-cluster index map
  all_labels <- c(S$atoms, Q$labels[1:Q$n_clusters])
  unique_labels <- unique(all_labels)
  label_to_cluster <- setNames(seq_along(unique_labels), unique_labels)
  
  # Step 6: Draw from posterior predictive distribution
  pred_samples <- numeric(N_predict)
  pred_clusters <- integer(N_predict)
  
  for (i in 1:N_predict) {
    if (runif(1) < R) {
      # From empirical component
      label <- sample(S$atoms, 1, prob = S$weights)
    } else {
      probs_Q <- Q$counts / Q$n_clusters
      label <- sample(Q$labels[1:Q$n_clusters], 1, prob = probs_Q)
    }
    
    pred_samples[i]  <- label
    pred_clusters[i] <- label_to_cluster[as.character(label)]
  }
  
  return(list(
    R              = R,
    S              = S,
    Q              = Q,
    pred_samples   = pred_samples,
    pred_clusters  = pred_clusters,
    label_map      = label_to_cluster
  ))
}

set.seed(123)
# Generate observations from a prior PY process
obs_data <- py_sample_lab(
  N            = 10,
  alpha        = 5,
  sigma        = 0.3,
  base_sampler = function() rnorm(1)
)

# Sample from posterior given those observations
post_res <- post_py_sample(
  N_predict = 100,
  alpha     = 5,
  sigma     = 0.3,
  data_py   = obs_data
)
