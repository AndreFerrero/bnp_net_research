py_pred = function(it, alpha, sigma){
  
  old_pi = function(x, this_x, alpha, sigma) {
    n = length(x)
    kn = length(unique(x))
    # keep track instead of recomputing
    nx = sum(this_x == x)
    
    old_x = (nx - sigma)/ (alpha + n)
    
    return(old_x)
  }
  
  new_pi = function(x, alpha, sigma) {
    n = length(x)
    kn = length(unique(x))
    
    new_x = (alpha + kn * sigma)/(alpha + n)
    
    return(new_x)
  }
  
  # initialise vector
  x = c(1)
  kn = 1
  
  t = 1
  
  while(t < it){
    old_p = c()
    p = c()
    
    # use unique to avoid repetitions
    # repetitions => more probabilities than needed in sample
    for(i in unique(x)){
      old_p = c(old_p, old_pi(x, i, alpha, sigma))
    }
    
    new_p = new_pi(x, alpha, sigma)
    
    p = c(old_p, new_p)
    
    # sample integer in the range of observed + new values
    new_x = sample(1:(kn+1), 1, replace = T, prob = p)
    
    # create static x before
    x = c(x, new_x)
    
    # update unique values
    kn = length(unique(x))
    
    t = t + 1
  }
  
  return(x)
}



py_net = function(N, alpha, sigma, plot = F){

  # t1 = Sys.time()
  x = py_pred(N, alpha[1], sigma[1])
  y = py_pred(N, alpha[2], sigma[2])
  # t2 = Sys.time()
  
  edges = data.frame(x = x,
                     y = y)
  tedges = table(edges) 
  
  # t3 = Sys.time()
  # 
  # print(t2-t1)
  # print(t3 - t2)
  
  return(list(adj_m = tedges,
              matr = matrix(tedges, ncol = ncol(tedges),
                               dimnames = dimnames(tedges))
              )
         )
}

py_net(1000, c(5,5), c(0,0))

