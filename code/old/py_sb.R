py_sb = function(it, alpha, sigma){
  w = c()
  wsum = 0
  csum = 0
  k = 1
  x = c()
  
  for(t in 1:it){
    u = runif(1)
    
    while(wsum < u){
      V = rbeta(1, 1 - sigma, alpha + k * sigma)
      
      wnew = (1 - wsum)*V
      
      w = c(w, wnew)
      wsum = wsum + wnew
      csum = c(csum, wsum)
      
      k = k + 1
    }
    
    thisx = max(which(csum < u))
    x = c(x, thisx)
  }
  
  return(list(w = w,
              x = x))
}