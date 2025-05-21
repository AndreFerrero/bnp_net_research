sim_net_Nsim = function(it, Nsim, m, alpha, sigma = c(0,0)){
  # m = motif number(s)
  # Nsim = Sample sizes to simulate
  
  library(bmotif)
  
  nets = list()
  motifs = list()
  m_freq = list()
  
  for(i in 1:length(Nsim)){
    nets_Nsim = list()
    motifs_Nsim = list()
    m_freq_Nsim = c()
    
    # for every Nsim draw it nets
    for(j in 1:it){
      nets_Nsim[[j]] = py_net(Nsim[i], alpha, sigma)$matr
      motifs_Nsim[[j]] = mcount(nets_Nsim[[j]], normalisation = F,
                                mean_weight = F, standard_dev = F)
      m_freq_Nsim = c(m_freq_Nsim,
                      subset(motifs_Nsim[[j]], motif == m)$frequency)
    }
    nets[[i]] = nets_Nsim
    motifs[[i]] = motifs_Nsim
    m_freq[[i]] = m_freq_Nsim
    
    
    print(paste("Nsim = ", Nsim[i]))}
  
  
  return(list(
    Nsim = Nsim,
    m_freq = m_freq
  ))
}

t1 = Sys.time()
sim3_Nsim_dp55 = sim_net_Nsim(100, seq(100, 1000, 100), 3, c(5,5))
t2 = Sys.time()
t2-t1

save(sim3_Nsim_dp55, file = "sim3_Nsim_dp55.Rdata")
boxplot(sim3_Nsim_dp55$m_freq, names = sim3_Nsim_dp55$Nsim,
        main = "DP, alpha = (5,5)")

sim3_Nsim_dp520 = sim_net_Nsim(100, seq(100, 1000, 100), 3, c(5,20))
save(sim3_Nsim_dp520, file = "sim3_Nsim_dp520.Rdata")
boxplot(sim3_Nsim_dp520$m_freq, names = sim3_Nsim_dp520$Nsim,
        main = "DP, alpha = (5,20)")
