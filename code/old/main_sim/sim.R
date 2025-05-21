sim = function(s, Nsim, m, alpha, sigma = c(0,0)){
  # m = motif number(s)
  # Nsim = Sample sizes to simulate
  
  library(bmotif)
  
  m_freq = list()
  time = numeric(length(Nsim))
  
  for(n in seq_along(Nsim)){
    
    t1 = Sys.time()
    m_freq[[n]] = replicate(s, net_n3(Nsim[n], m, alpha))
    t2 = Sys.time()
    
    time[n] = difftime(t2, t1, units = "secs") |> as.numeric()
    
    print(paste("Nsim = ", Nsim[n]))}
  
  
  return(list(m_freq = m_freq,
              time = time))
}





# sizes = seq(500, 10000, 500)
# 
# sim3510 = sim_net_Nsim(s = 100,
#                        sizes,
#                        m = c(3,5,10),
#                        alpha = c(5,5))
# 
# save(sim3510, file = "sims/sim3510.Rdata")


# sizes = seq(500, 100000, 500)
# 
# sim_to_100000 = sim(s = 50,
#                              sizes,
#                              m = c(3,5,10),
#                              alpha = c(5,5))
# save(sim_to_100000, file = "sims/sim_to_100000.Rdata")

# t1 = Sys.time()
# # net_n3(100000, c(3, 5), c(5,5))
# replicate(100, net_n3(100000, c(3, 5), c(5,5)))
# t2 = Sys.time()
# 
# difftime(t2,t1, units = "min")

