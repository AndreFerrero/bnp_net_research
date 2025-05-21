
m_cond_dist_parallel <- function(N, t, d, alpha) {
  # Detect the number of cores available and register them for parallel execution
  num_cores <- detectCores() - 1  # Leave one core free for system tasks
  registerDoParallel(cores = num_cores)
  
  # Parallelize the simulations using foreach
  results <- foreach(i = seq_len(t),
                     .packages = c("ggplot2", "dplyr"),
                     .export = c("py_sample", "py_predictive_draw3", "m_xi"),
                     .combine = rbind) %dopar% {
                       xA <- py_sample(N, alpha)
                       xB <- py_sample(N, alpha)
                       
                       m_values <- numeric(d)
                       for (j in 1:d) {
                         xiA <- py_pred(alpha, 0, length(xA$x), xA$counts_x, xA$idx_list_x)
                         xiB <- py_pred(alpha, 0, length(xB$x), xB$counts_x, xB$idx_list_x)
                         m_values[j] <- m_xi(c(xiA, xiB), xA, xB)
                       }
                       
                       # Return a properly formatted data frame
                       data.frame(
                         Simulation = factor(rep(i, d)),  # Repeat simulation identifier for each value
                         Value = m_values  # Numeric values
                       )
                     }
  
  # Stop the parallel backend once done
  stopImplicitCluster()
  
  return(results)  # This is now a single clean dataframe
}
