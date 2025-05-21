net_n3 = function(N, m, alpha, sigma = c(0, 0)) {
  
  alpha1 = alpha[1]
  alpha2 = alpha[2]
  
  sigma1 = sigma[1]
  sigma2 = sigma[2]
  
  x = y = numeric(N)
  x[1] = y[1] = 1
  
  counts_x = counts_y = numeric(N)
  counts_x[1] = counts_y[1] = 1
  
  # Initialize idx_list for active clusters (non-zero entries)
  idx_list_x = idx_list_y = c(1)  # Only the first cluster (index 1) is initially non-zero
  
  for (i in 2:N) {
    # Sampling for x and y, passing the active cluster indices (idx_list_x, idx_list_y)
    x[i] <- py_pred(alpha1, sigma1, i, counts_x, idx_list_x)
    y[i] <- py_pred(alpha2, sigma2, i, counts_y, idx_list_y)
    
    # Update counts and active indices for x
    counts_x[x[i]] = counts_x[x[i]] + 1
    if (counts_x[x[i]] == 1) {  # If this was the first time this cluster appeared
      idx_list_x = c(idx_list_x, x[i])  # Add the new active cluster to the list
    }
    
    # Update counts and active indices for y
    counts_y[y[i]] = counts_y[y[i]] + 1
    if (counts_y[y[i]] == 1) {  # If this was the first time this cluster appeared
      idx_list_y = c(idx_list_y, y[i])  # Add the new active cluster to the list
    }
  }
  
  inter <- cbind(x, y)
  network <- as.matrix(table(inter[, 1], inter[, 2]))
  network <- matrix(network, ncol = ncol(network))
  
  freqs = mcount(network, normalisation = F,
                 mean_weight = F, standard_dev = F)[, "frequency"]
  
  return(freqs[m])
}

net_n3(100000, c(3, 5), c(5,5))