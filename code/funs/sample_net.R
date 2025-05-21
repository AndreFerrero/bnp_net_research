sample_net = function(N, alpha, sigma = c(0, 0)) {
  # N = number of edges
  
  alphaA = alpha[1]
  alphaB = alpha[2]
  
  sigmaA = sigma[1]
  sigmaB = sigma[2]
  
  xA = py_sample(N, alphaA, sigmaA)
  xB = py_sample(N, alphaB, sigmaB)
  
  edges = data.frame(X_A = xA$x,
                     X_B = xB$x)
  
  result = list(edges = edges,
                xA = xA,
                xB = xB)
  
  return(result)
}
