py_sample_lab <- function(N, alpha = 5, sigma = 0,
                                  base_sampler = function() rnorm(1)) {
  # x[i]           : which cluster observation i is in
  # labels_list[k] : the label (e.g. a real number) associated to cluster k
  # label_of_obs[i]: the label drawn from labels_list corresponding to x[i]
  
  x           <- integer(N)
  label_of_obs<- numeric(N)
  
  # initialize with one cluster
  x[1]            <- 1
  counts_x        <- integer(N); counts_x[1] <- 1
  idx_list_x      <- 1     # active cluster indices
  labels_list     <- base_sampler()  # one label for the first cluster
  label_of_obs[1] <- labels_list[1]
  
  for (i in 2:N) {
    k <- length(idx_list_x)
    # compute predictive probs for existing clusters...
    probs <- c(counts_x[idx_list_x] - sigma,
               alpha + sigma * k) / (alpha + i - 1)
    
    # sample a cluster index in 1:(k+1)
    choice <- sample.int(k + 1, 1, prob = probs)
    
    if (choice <= k) {
      # existing cluster
      clust <- idx_list_x[choice]
    } else {
      # new cluster: assign next integer k+1
      clust <- k + 1
      idx_list_x <- c(idx_list_x, clust)
      # draw a fresh label from the base (e.g. N(0,1)) and store it
      labels_list[clust] <- base_sampler()
    }
    
    # record assignment and update counts
    x[i]            <- clust
    counts_x[clust] <- counts_x[clust] + 1
    
    # set the observed‐label for this data point
    label_of_obs[i] <- labels_list[clust]
  }
  
  list(
    clusters     = x,
    counts       = counts_x[idx_list_x],
    labels       = labels_list,
    obs_labels   = label_of_obs,
    n_clusters   = length(idx_list_x)
  )
}

## Example usage:
# set.seed(42)
# res <- py_sample_lab(
#   N          = 10,
#   alpha      = 2,
#   sigma      = 0.5,
#   base_sampler = function() rnorm(1, mean = 0, sd = 1)
# )
# 
# str(res)
