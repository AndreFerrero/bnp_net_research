sample_net = function(N, alpha, sigma = c(0, 0)) {
  # N = number of edges
  
  alpha1 = alpha[1]
  alpha2 = alpha[2]
  
  sigma1 = sigma[1]
  sigma2 = sigma[2]
  
  x = y = numeric(N)
  x[1] = y[1] = 1
  
  counts_x = counts_y = numeric(N)
  counts_x[1] = counts_y[1] = 1
  
  idx_list_x = idx_list_y = c(1)  # active cluster indices
  
  for (i in 2:N) {
    # Sample cluster index for node group X and Y
    x[i] <- py_pred(alpha1, sigma1, i, counts_x, idx_list_x)
    y[i] <- py_pred(alpha2, sigma2, i, counts_y, idx_list_y)
    
    # Update counts and idx list for X
    counts_x[x[i]] = counts_x[x[i]] + 1
    if (!(x[i] %in% idx_list_x)) {
      idx_list_x = c(idx_list_x, x[i])
    }
    
    # Update counts and idx list for Y
    counts_y[y[i]] = counts_y[y[i]] + 1
    if (!(y[i] %in% idx_list_y)) {
      idx_list_y = c(idx_list_y, y[i])
    }
  }
  
  # Create edge list
  edges <- data.frame(X_A = x, X_B = y)
  
  # Prepare detailed output similar to py_sample_lab
  result <- list(
    edges = edges,
    X = list(
      clusters = x,
      counts = counts_x[idx_list_x],
      cluster_ids = idx_list_x,
      n_clusters = length(idx_list_x)
    ),
    Y = list(
      clusters = y,
      counts = counts_y[idx_list_y],
      cluster_ids = idx_list_y,
      n_clusters = length(idx_list_y)
    )
  )
  
  return(result)
}
