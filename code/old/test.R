library(parallel)

sim_net_Nsim_parallel <- function(it, Nsim, motifs, alpha,
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
    results = parLapply(cl, 1:it, function(j){
      net = py_net(Nsim[i], alpha, sigma)$matr
      motif_data = mcount(net, normalisation = FALSE,
                          mean_weight = FALSE, standard_dev = FALSE)
      list(net = net, motif_data = motif_data)
    })
    
    motifs_matr = lapply(results, `[[`, "motif_data")
    
    # motif_freq = lapply(motifs_matr, function {subset(motifs_matr, motif == m)$frequency})
    # 
    # for(m in seq_along(motifs)){
    #   motif_freq = 
    # }
    
    m_freq[[i]][[m]] = unlist(lapply(results, `[[`, "motif_freq"))
  }
    
  print(paste("Nsim =", Nsim[i]))
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  return(list(
    Nsim = Nsim,
    m_freq = m_freq
  ))
}

s = sim_net_Nsim_parallel(4, 100:105, c(3, 5), c(5,5))

