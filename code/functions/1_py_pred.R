py_pred <- function(alpha, sigma, n, counts, idx_list) {
  k <- length(idx_list)  
  prob <- c(counts[idx_list] - sigma, alpha + sigma * k) / (alpha + n)
  pred <- sample(c(idx_list, k + 1), 1, prob = prob)
  return(pred)
}