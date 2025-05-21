library(parallel)

sim <- function(it, Nsim, m, alpha,
                                  sigma = c(0,0),
                                  ncores = detectCores() - 1) {
  library(bmotif)
  # m = motif number(s)
  # Nsim = Sample sizes to simulate
  # ncores = Number of CPU cores to use (default: all except one)
  
  # Create a cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("py_net", "py_pred", "mcount"))
  
  m_freq = vector("list", length(Nsim))
  
  for(i in seq_along(Nsim)) {
    # Parallelized computation using parLapply
    results = parLapply(cl, 1:it, function(j) {
      net = py_net(Nsim[i], alpha, sigma)$matr
      motif_data = mcount(net, normalisation = FALSE,
                          mean_weight = FALSE, standard_dev = FALSE)
      motif_freq = subset(motif_data, motif == m)$frequency
      list(net = net, motif_data = motif_data, motif_freq = motif_freq)
    })
    
    m_freq[[i]] = unlist(lapply(results, `[[`, "motif_freq"))
    
    print(paste("Nsim =", Nsim[i]))
  }
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  return(list(
    Nsim = Nsim,
    m_freq = m_freq
  ))
}

sim_net_Nsim_parallel(4, 100:105, 3, c(5,5))
