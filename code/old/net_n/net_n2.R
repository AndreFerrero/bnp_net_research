py_predictive_draw2 <- function(alpha, sigma, n, counts) {
  counts = counts[counts > 0]
  k = length(counts)
  
  prob <- c(counts - sigma, alpha + sigma * k)/(alpha + n)
  pred <- sample(1:(k+1), 1, prob=prob)
  return(pred)
}

net_n2 = function(N, m, alpha, sigma = c(0,0)){
  
  alpha1 = alpha[1]
  alpha2 = alpha[2]
  
  sigma1 = sigma[1]
  sigma2 = sigma[2]
  
  x = y = numeric(N)
  x[1] = y[1] = 1
  n = 1
  
  counts_x = counts_y = numeric(N)
  counts_x[1] = counts_y[1] = 1
  
  for (i in 2:N) {
    x[i] <- py_predictive_draw2(alpha1, sigma1, n, counts_x)
    y[i] <- py_predictive_draw2(alpha2, sigma2, n, counts_y)
    n = n + 1
    counts_x[x[i]] = counts_x[x[i]] + 1
    counts_y[y[i]] = counts_y[y[i]] + 1
  }
  
  inter <- cbind(x,y)
  network <- as.matrix(table(inter[,1], inter[,2]))
  network <- matrix(network, ncol = ncol(network))
  
  freqs = mcount(network, normalisation = F,
                 mean_weight = F, standard_dev = F)[,"frequency"]
  
  return(freqs[m])
}

net_n2(100000, c(3, 5), c(5,5))
