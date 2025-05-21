sim_net = function(Nsim, m, alpha, sigma = c(0,0)){
  # m = motif number(s)
  # Nsim = Sample sizes to simulate
  
  library(bmotif)
  
  nets = list()
  motifs = list()
  m_freq = c()
  
  for(i in 1:length(Nsim)){
    nets[[i]] = py_net(Nsim[i], alpha, sigma)$matr
    motifs[[i]] = mcount(nets[[i]], normalisation = F,
                         mean_weight = F, standard_dev = F)
    
    m_freq = c(m_freq,
               subset(motifs[[i]], motif == m)$frequency)
    
    print(paste("Nsim = ", Nsim[i]))
  }
  
  
  return(data.frame(
    Nsim = Nsim,
    m_freq = m_freq
  ))
}

