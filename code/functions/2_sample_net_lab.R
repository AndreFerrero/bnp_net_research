sample_net_lab <- function(N,
                           alpha = c(5, 5),
                           sigma = c(0, 0),
                           base_sampler = function() rnorm(1)) {
  # N           : number of edges to generate
  # alpha       : vector length 2; concentration for A- and B-side CRPs
  # sigma       : vector length 2; discount for A- and B-side CRPs
  # base_sampler: function returning one draw from the base measure
  
  if (length(alpha) != 2 || length(sigma) != 2) {
    stop("alpha and sigma must be numeric vectors of length 2")
  }
  
  # Sample clusters & labels for A-side
  g_A <- py_sample_lab(
    N           = N,
    alpha       = alpha[1],
    sigma       = sigma[1],
    base_sampler= base_sampler
  )
  
  # Sample clusters & labels for B-side
  g_B <- py_sample_lab(
    N           = N,
    alpha       = alpha[2],
    sigma       = sigma[2],
    base_sampler= base_sampler
  )
  
  # Build the bipartite edge list
  edges <- data.frame(
    X_A = g_A$clusters,
    X_B = g_B$clusters
  )
  
  # Return edges plus the full CRP outputs for each side
  list(
    edges = edges,
    X = g_A,
    Y = g_B
  )
}
