plot_sim = function(sim, m, alpha, extra = F){
  
  plot(sim$Nsim, sim$m_freq, pch = 19,
       main = paste0("Simulation for motif", m,
                     " alpha = ", "(", alpha[1], ',', alpha[2], ')'))
  
  if(extra){
    plot(sim$Nsim, log(sim$m_freq), pch = 19,
         main = paste0("Simulation for motif", m, " (log scale)",
                       " alpha = ", "(", alpha[1], ',', alpha[2], ')'))
    
    plot(log(sim$Nsim), log(sim$m_freq), pch = 19,
         main = paste0("Simulation for motif", m, " (log-log scale)",
                       " alpha = ", "(", alpha[1], ',', alpha[2], ')'))
  }
}
