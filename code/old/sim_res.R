
## Motif 3 ####
### 5 5 ####
p_sim3_Nsim_dp55_to4000 = sim_net_Nsim_parallel(50,
                                                seq(100, 4000, 200),
                                                3,
                                                c(5,5))
save(p_sim3_Nsim_dp55_to4000, file = "p_sim3_Nsim_dp55_to4000.Rdata")

boxplot(p_sim3_Nsim_dp55_to4000$m_freq, names = p_sim3_Nsim_dp55_to4000$Nsim,
        main = "DP, alpha = (5,5), Motif 3")
means <- sapply(p_sim3_Nsim_dp55_to4000$m_freq, mean)
points(1:length(means), means, col = "red", pch = 19)

plot(log(p_sim3_Nsim_dp55_to4000$Nsim),
     means, pch = 19)



### 5 20 ####
p_sim3_Nsim_dp520_to4000 = sim_net_Nsim_parallel(50,
                                                 seq(100, 4000, 200),
                                                 3,
                                                 c(5,20))
save(p_sim3_Nsim_dp520_to4000, file = "p_sim3_Nsim_dp520_to4000.Rdata")

boxplot(p_sim3_Nsim_dp520_to4000$m_freq, names = p_sim3_Nsim_dp520_to4000$Nsim,
        main = "DP, alpha = (5,20), Motif 3")
means <- sapply(p_sim3_Nsim_dp520_to4000$m_freq, mean)
points(1:length(means), means, col = "red", pch = 19)

plot(log(p_sim3_Nsim_dp520_to4000$Nsim),
     means, pch = 19)


## Motif 5 ####
p_sim5_Nsim_dp55 = sim_net_Nsim_parallel(100, seq(100, 1000, 100), 5, c(5,5))
save(p_sim5_Nsim_dp55, file = "p_sim5_Nsim_dp55.Rdata")

boxplot(p_sim5_Nsim_dp55$m_freq, names = p_sim5_Nsim_dp55$Nsim,
        main = "DP, alpha = (5,5), Motif 5")
means <- sapply(p_sim5_Nsim_dp55$m_freq, mean)
points(1:length(means), means, col = "red", pch = 19)




## Motif 10 ####
p_sim10_Nsim_dp55 = sim_net_Nsim_parallel(100, seq(100, 1000, 100), 10, c(5,5))
save(p_sim10_Nsim_dp55, file = "p_sim10_Nsim_dp55.Rdata")

boxplot(p_sim10_Nsim_dp55$m_freq, names = p_sim10_Nsim_dp55$Nsim,
        main = "DP, alpha = (5,5), Motif 10")
means <- sapply(p_sim10_Nsim_dp55$m_freq, mean)
points(1:length(means), means, col = "red", pch = 19)


