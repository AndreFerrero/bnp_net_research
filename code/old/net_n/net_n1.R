py_predictive_draw1 <- function(alpha, sigma, given, n) {
  counts = tabulate(given)
  k = length(counts)
  
  prob <- c(counts - sigma, alpha + sigma * k)/(alpha + n)
  pred <- sample(1:(k+1), 1, prob=prob)
  return(pred)
}

net_n1 = function(N, m, alpha, sigma = c(0,0)){

  alpha1 = alpha[1]
  alpha2 = alpha[2]

  sigma1 = sigma[1]
  sigma2 = sigma[2]

  x <- rep(0,N)
  y <- rep(0,N)
  x[1] <- 1
  y[1] <- 1
  n = 1

  for (i in 2:N) {
    x[i] <- py_predictive_draw1(alpha1, sigma1, x[1:(i-1)], n)
    y[i] <- py_predictive_draw1(alpha2, sigma2, y[1:(i-1)], n)
    n = n + 1
  }

  inter <- cbind(x,y)
  network <- as.matrix(table(inter[,1], inter[,2]))
  network <- matrix(network, ncol = ncol(network))

  freqs = mcount(network, normalisation = F,
                 mean_weight = F, standard_dev = F)[,"frequency"]

  return(freqs[m])
}

# net_n(100000, c(3, 5), c(5,5))
