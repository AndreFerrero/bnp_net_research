py_sample <- function(N, alpha = 5, sigma = 0) {
  x <- numeric(N)
  x[1] <- 1
  
  counts_x <- numeric(N)
  counts_x[1] <- 1
  
  idx_list_x <- c(1)
  
  for (i in 2:N) {
    k <- length(idx_list)  
    prob <- c(counts[idx_list] - sigma, alpha + sigma * k) / (alpha + n)
    x[i] <- sample(c(idx_list, k + 1), 1, prob = prob)
    counts_x[x[i]] <- counts_x[x[i]] + 1
    if (counts_x[x[i]] == 1) {  
      idx_list_x <- c(idx_list_x, x[i])
    }
  }
  
  return(list(x = x, counts_x = counts_x, idx_list_x = idx_list_x,
              K = length(idx_list_x)))
}