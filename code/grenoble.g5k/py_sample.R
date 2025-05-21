py_sample <- function(N, alpha = 5, sigma = 0) {
  x <- numeric(N)
  x[1] <- 1
  
  counts <- numeric(N)
  counts[1] <- 1
  
  active_nodes <- c(1)
  
  for (n in 2:N) {
    k <- length(active_nodes)  
    prob <- c(counts[active_nodes] - sigma, alpha + sigma * k)
    x[n] <- sample(c(active_nodes, k + 1), 1, prob = prob)
    counts[x[n]] <- counts[x[n]] + 1
    if (counts[x[n]] == 1) {  
      active_nodes <- c(active_nodes, x[n])
    }
  }
  
  return(list(x = x,
              counts = counts,
              active_nodes = active_nodes,
              active_counts = counts[active_nodes],
              K = length(active_nodes)))
}