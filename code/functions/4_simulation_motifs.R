#############################
# 1. Sequential simulation #
#############################

simulation_networks_loop <- function(S, N, 
                                     alphaA = 5, sigmaA = 0, 
                                     alphaB = 5, sigmaB = 0) {
  # Initialize matrices to store motif counts, indicator values and cumulative motif counts.
  M_mat <- matrix(0, nrow = S, ncol = N)
  I_mat <- matrix(0, nrow = S, ncol = N)
  cumulativeMotifs_mat <- matrix(0, nrow = S, ncol = N)
  
  X_A <- vector("list", S)
  X_B <- vector("list", S)
  
  cat("Starting sequential simulations...\n")
  
  for (s in 1:S) {
    if (s %% (S / 10) == 0) {  # Print progress every 10% completion
      cat(sprintf("Progress: %d%% completed (%d/%d simulations)\n", 
                  round(100 * s / S), s, S))
    }
    
    # Generate sequences for group A and group B using the PY/DP process.
    resA <- py_sample(N, alphaA, sigmaA)
    resB <- py_sample(N, alphaB, sigmaB)
    xA <- resA$x
    xB <- resB$x
    
    # Count motifs using the count_motifs function.
    motif_counts <- count_motifs(xA, xB)
    M_mat[s, ] <- motif_counts$M
    I_mat[s, ] <- motif_counts$Indicator
    cumulativeMotifs_mat[s, ] <- cumsum(motif_counts$M)
    
    X_A[[s]] <- xA
    X_B[[s]] <- xB
  }
  
  cat("Sequential simulations completed.\n")
  
  # Compute averages over simulations for each network size (draw)
  avgM <- colMeans(M_mat)                     # Average motif count per draw
  avgZero <- colMeans(I_mat)                    # Average indicator (1 means no motif added)
  avgCumulativeMotifs <- colMeans(cumulativeMotifs_mat)  # Average cumulative motifs
  
  return(list(M = M_mat, Indicator = I_mat, 
              cumulativeMotifs = cumulativeMotifs_mat,
              X_A = X_A, X_B = X_B,
              avgM = avgM, avgZero = avgZero,
              avgCumulativeMotifs = avgCumulativeMotifs))
}


##########################################
# 2. Parallelized simulation over S (S)  #
##########################################

simulation_networks_parallel <- function(S, N, 
                                         alphaA = 5, sigmaA = 0, 
                                         alphaB = 5, sigmaB = 0, 
                                         ncores = 2) {
  require(parallel)
  
  cl <- makeCluster(ncores)
  clusterExport(cl, c("py_sample", "count_motifs", "py_predictive_draw3"))
  
  cat("Starting parallel simulations using", ncores, "cores...\n")
  
  one_simulation <- function(sim_id) {
    cat(sprintf("Core %d: Running simulation %d/%d...\n", Sys.getpid(), sim_id, S))
    resA <- py_sample(N, alphaA, sigmaA)
    resB <- py_sample(N, alphaB, sigmaB)
    
    xA <- resA$x
    xB <- resB$x
    motif_counts <- count_motifs(xA, xB)
    cumulativeMotifs <- cumsum(motif_counts$M)
    
    return(list(M = motif_counts$M, 
                Indicator = motif_counts$Indicator, 
                cumulativeMotifs = cumulativeMotifs,
                X_A = xA, X_B = xB))
  }
  
  sim_list <- parLapply(cl, 1:S, one_simulation)
  stopCluster(cl)
  
  cat("Parallel simulations completed.\n")
  
  M_mat <- t(sapply(sim_list, function(res) res$M))
  I_mat <- t(sapply(sim_list, function(res) res$Indicator))
  cumulativeMotifs_mat <- t(sapply(sim_list, function(res) res$cumulativeMotifs))
  X_A <- lapply(sim_list, function(res) res$X_A)
  X_B <- lapply(sim_list, function(res) res$X_B)
  
  avgM <- colMeans(M_mat)
  avgZero <- colMeans(I_mat)
  avgCumulativeMotifs <- colMeans(cumulativeMotifs_mat)
  
  return(list(M = M_mat, Indicator = I_mat, 
              cumulativeMotifs = cumulativeMotifs_mat,
              X_A = X_A, X_B = X_B,
              avgM = avgM, avgZero = avgZero,
              avgCumulativeMotifs = avgCumulativeMotifs))
}