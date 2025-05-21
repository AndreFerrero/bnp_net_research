sim_par2 <- function(s, Nsim, m, alpha, sigma = c(0,0)){
  # Load necessary libraries
  library(bmotif)
  library(parallel)
  
  # Detect number of available cores
  n_cores <- detectCores() - 1 # Use one less than max cores to avoid system overload
  
  # Function to run a single simulation
  sim_fun <- function(n) {
    t1 <- Sys.time()
    m_freq <- replicate(s, net_n3(n, m, alpha))
    t2 <- Sys.time()
    time_taken <- as.numeric(difftime(t2, t1, units = "secs"))
    print(paste("Nsim =", n))
    return(list(m_freq = m_freq, time = time_taken))
  }
  
  # Initialize cluster
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("py_predictive_draw3", "net_n3", "mcount"))
  
  # Run simulations in parallel
  results <- parLapply(cl, Nsim, sim_fun)
  
  # Stop cluster
  stopCluster(cl)
  
  # Extract results
  m_freq <- lapply(results, `[[`, "m_freq")
  time <- sapply(results, `[[`, "time")
  
  return(list(m_freq = m_freq, time = time))
}

sizes = seq(500, 100000, 500)

sim_to_100000 = sim_par2(s = 50,
                             sizes,
                             m = c(3,5,10),
                             alpha = c(5,5))
save(sim_to_100000, file = "sims/sim_to_100000.Rdata")