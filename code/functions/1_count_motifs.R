# This function counts motifs by comparing the current edge with previous unique connections.
# Note: This function counts only unique A nodes connected to the same B node.
count_motifs <- function(xA, xB) {
  N <- length(xA)  
  M <- numeric(N)        # Number of new motifs added at each draw
  Indicator <- numeric(N)  # 1 if no new motif, 0 if a motif was added
  
  # For the first draw: by definition no motif is created.
  M[1] <- 0
  Indicator[1] <- 1  
  
  observed_edges <- list()
  observed_edges[[1]] <- paste(xA[1], xB[1], sep = "-")
  
  for (n in 2:N) {
    xiA <- xA[n]
    xiB <- xB[n]
    current_edge <- paste(xiA, xiB, sep = "-")
    
    if (!(current_edge %in% observed_edges)) {
      # Extract the subset of previous xA values where the corresponding xB equals xiB
      previous_As <- xA[1:(n-1)][xB[1:(n-1)] == xiB]
      unique_As <- unique(previous_As)
      count_unique_As <- length(unique_As)
      
      if (count_unique_As > 0) {
        M[n] <- count_unique_As
        Indicator[n] <- 0  # 0 indicates that a motif was added at this draw
      } else {
        Indicator[n] <- 1  # 1 indicates no motif was added
      }
      
      # Add this edge to the observed list so it isn't counted again
      observed_edges[[length(observed_edges) + 1]] <- current_edge
    } else {
      Indicator[n] <- 1
    }
  }
  
  return(list(M = M, Indicator = Indicator))
}