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
  clusterExport(cl, c("net_n"))
  
  m_freq = vector("list", length(Nsim))
  
  for(n in seq_along(Nsim)) {
    results = parLapply(cl, 1:it, net_n(n,m,alpha))
    
    
    
    print(paste("Nsim =", Nsim[n]))
  }
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  return(results)
}

sim(4, 100:105, c(3,5), c(5,5))
