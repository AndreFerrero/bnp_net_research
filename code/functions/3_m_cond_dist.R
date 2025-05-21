# Function to run the simulations and return a dataframe
m_cond_dist <- function(N, t, d, alpha) {
  results <- data.frame()
  
  for (i in seq_len(t)) {
    xA <- py_sample(N, alpha)
    xB <- py_sample(N, alpha)
    
    m_values <- replicate(d, {
      xiA <- py_pred(alpha, 0, length(xA$x), xA$counts_x, xA$idx_list_x)
      xiB <- py_pred(alpha, 0, length(xB$x), xB$counts_x, xB$idx_list_x)
      m_xi(c(xiA, xiB), xA, xB)
    }, simplify = TRUE)  # Collect values in a vector
    
    # Convert to dataframe
    temp_df <- data.frame(
      Simulation = factor(i),  # Simulation identifier (1, 2, 3, ...)
      Value = m_values  # Values for histogram
    )
    
    results <- bind_rows(results, temp_df)  # Combine all results
  }
  
  return(results)
}