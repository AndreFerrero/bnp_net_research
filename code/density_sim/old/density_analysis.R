
# ----------------------------
# 1. Generate the networks
# ----------------------------
# set.seed(123)  # for reproducibility in the master process
# nsim <- 100      # number of simulations
# N <- 1000000      # number of edges per simulation
# 
# # Setup a cluster using available cores minus one
# numCores <- detectCores() - 1
# cl <- makeCluster(numCores)
# 
# # Export the necessary function and variables to the workers.
# # Ensure sample_net is defined/available in this session.
# clusterExport(cl, varlist = c("py_predictive_draw3",
#                               "sample_net", "N"))
# 
# # Run simulations in parallel using parLapply.
# # Each worker generates one network returning a dataframe with columns X_A and X_B.
# sim_edgesDP2020 <- parLapply(cl, 1:nsim, function(i) {
#   # Optionally, use set.seed(i) here if you need reproducibility across workers.
#   sample_net(N, alpha = c(20,20), sigma = c(0,0))
# })
# 
# # Stop the cluster to free resources
# stopCluster(cl)

save(sim_edgesDP2020,
     file = "C:/Users/anferrer/OneDrive - Alma Mater Studiorum Università di Bologna/University/Internship/research/sims/dp2020.Rdata")

# ----------------------------
# 2. Static Analysis: Summary of Each Full Network
# ----------------------------
# For each simulation, compute:
# - The number of unique edges (removing any duplicates, if present)
# - The number of unique nodes in each partition
# - The maximum possible edges (nA * nB)
# - The density = unique_edges / (nA * nB)
# 
sim_edges = sim_edgesPY552525
sim_edges = sim_edgesPY557575

# 
# static_results <- data.frame(Simulation = integer(0),
#                              UniqueEdges = integer(0),
#                              nA = integer(0),
#                              nB = integer(0),
#                              MaxPossible = integer(0),
#                              Density = numeric(0))
# 
# nsim <- length(sim_edges)
# 
# for (i in 1:nsim) {
#   full_net <- sim_edges[[i]]
#   
#   # Consider only unique edges since your model might allow repeats
#   unique_net <- unique(full_net)
#   n_unique <- nrow(unique_net)
#   
#   # Unique nodes in each partition:
#   nA <- length(unique(unique_net$X_A))
#   nB <- length(unique(unique_net$X_B))
#   
#   max_possible <- nA * nB
#   density <- n_unique / max_possible
#   
#   static_results <- rbind(static_results,
#                           data.frame(Simulation = i,
#                                      UniqueEdges = n_unique,
#                                      nA = nA,
#                                      nB = nB,
#                                      MaxPossible = max_possible,
#                                      Density = density))
# }
# 
# # Print the static summary results
# print(static_results)

# ----------------------------
# 3. Dynamic Analysis: Density Evolution with Thresholds
# ----------------------------
# We will "cut" the network at a set of predetermined thresholds (i.e. considering the first T edges)
thresholds <- c(100:500, seq(1000,100000, by = 1000))
nsim <- length(sim_edges)

# Create a list to store results for the dynamic analysis.
dyn_results_list <- list()

for (i in 1:nsim) {
  # Get the full network from the simulation list
  full_net <- sim_edges[[i]]
  
  # Loop over each threshold
  for (th in thresholds) {
    # Take the first 'th' edges
    sub_net <- full_net[1:th, ]
    
    # Ensure uniqueness
    sub_unique <- unique(sub_net)
    n_sub_unique <- nrow(sub_unique)
    
    # Count unique nodes in each partition from the subnetwork
    nA_sub <- length(unique(sub_unique$X_A))
    nB_sub <- length(unique(sub_unique$X_B))
    
    # Maximum possible edges for this subnetwork
    max_possible_sub <- nA_sub * nB_sub
    
    # Calculate density for this threshold
    density_sub <- n_sub_unique / max_possible_sub
    
    # Save the results into a data frame
    dyn_results_list[[length(dyn_results_list) + 1]] <- data.frame(
      Simulation = i,
      Threshold = th,
      UniqueEdges = n_sub_unique,
      nA = nA_sub,
      nB = nB_sub,
      MaxPossible = max_possible_sub,
      Density = density_sub
    )
  }
}

# Combine dynamic results into one data frame
dyn_results <- do.call(rbind, dyn_results_list)

# dens_PY75 <- dyn_results
# dens_folder <- here("res", "density_PY")
# save(dens_PY75, file = here(dens_folder, "dens_PY_0.75.Rdata"))

dens_PY25 <- dyn_results
dens_folder <- here("res", "density_PY")
save(dens_PY25, file = here(dens_folder, "dens_PY_0.25.Rdata"))

# ----------------------------
# 4. Plotting: Boxplots of Density Across Thresholds
# ----------------------------
# Plot a boxplot to show the distribution of density across simulations for each threshold

# ggplot(dyn_results, aes(x = factor(Threshold), y = Density)) +
#   geom_boxplot(fill = "skyblue", color = "darkblue", outlier.color = "red") +
#   labs(x = "Network Size",
#        y = "Density",
#        title = "DP(5) densities, 100 networks") +
#   theme_minimal(base_size = 14)

# ----------------------------
# 5. Plotting: Jittered Points for Each Threshold

# ggplot(dyn_results, aes(x = factor(Threshold), y = Density)) +
#   geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA, alpha = 0.5) +
#   geom_jitter(width = 0.3, height = 0, size = 1.5, alpha = 0.6, color = "black") +
#   labs(x = "Network Size",
#        y = "Density",
#        title = "DP(5) densities, 100 networks") +
#   theme_minimal(base_size = 14)
# ggsave("sims/pics/density_dp0505.png", bg = "white")
# 
# ggplot(dyn_results, aes(x = factor(Threshold), y = UniqueEdges)) +
#   geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA, alpha = 0.5) +
#   geom_jitter(width = 0.3, height = 0, size = 1.5, alpha = 0.6, color = "black") +
#   labs(x = "Network Size",
#        y = "Unique edges",
#        title = "PY(5, 0.6), 20 networks") +
#   theme_minimal(base_size = 14)
# ggsave("sims/pics/En_py5566.png", bg = "white")
# 
# ggplot(dyn_results, aes(x = factor(Threshold), y = MaxPossible)) +
#   geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA, alpha = 0.5) +
#   geom_jitter(width = 0.3, height = 0, size = 1.5, alpha = 0.6, color = "black") +
#   labs(x = "Network Size",
#        y = "NA x NB",
#        title = "PY(5, 0.6), 20 networks") +
#   theme_minimal(base_size = 14)
# ggsave("sims/pics/NAxNB_py5566.png", bg = "white")

ggplot(dens_PY75, aes(x = log(Threshold), y = log(Density),
                        group = factor(Simulation),
                        color = factor(Simulation))) +
  geom_line(alpha = 0.4, size = 0.2) +
  scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
  labs(x = expression(log(e(Y[n]))),
       y = expression(log(d(H[n]))),
       color = "Simulation") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

library(here)
dens_pics <- here("res", "pics", "density_analysis")
ggsave(filename = here(dens_pics, "PY_5__0.25__density.png"),
       width = 11,
       height = 8,
       dpi = 600,
       bg = "white")

# ggplot(dyn_results, aes(x = log(Threshold), y = log(UniqueEdges),
#                         group = factor(Simulation),
#                         color = factor(Simulation))) +
#   geom_line(alpha = 0.3, size = 0.8) +
#   scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
#   labs(x = "Observed Edges",
#        y = "Unique edges",
#        color = "Simulation") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none")
# 
# ggplot(dyn_results, aes(x = log(Threshold), y = log(MaxPossible),
#                         group = factor(Simulation),
#                         color = factor(Simulation))) +
#   geom_line(alpha = 0.3, size = 0.8) +
#   scale_color_viridis_d(option = "plasma", begin = 0, end = 1) +
#   labs(x = "Observed Edges",
#        y = "NAxNB",
#        color = "Simulation") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "none")
