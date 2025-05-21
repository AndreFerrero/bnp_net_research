library(here)

dens_folder <- here("dens")

load(here(dens_folder, "py_552525.Rdata"))
load(here(dens_folder, "py_5555.Rdata"))
load(here(dens_folder, "py_557575.Rdata"))

dens <- function(edges, ncores = parallel::detectCores(logical = FALSE)) {
  require(parallel)
  
  N    <- seq(50, 100000, by = 50)
  nsim <- length(edges)
  
  message("Launching density analysis with ", nsim,
          " simulations across ", ncores, " cores.")
  
  # mclapply will fork worker processes on Unix-alikes.
  results_list <- mclapply(
    X         = seq_len(nsim),
    FUN       = function(i) {
      message(sprintf("[sim %d/%d] START", i, nsim))
      
      full_net <- edges[[i]]
      sim_res  <- vector("list", length(N))
      
      for (j in seq_along(N)) {
        n <- N[j]
        if (n > nrow(full_net)) next
        
        sub_unique   <- unique(full_net[1:n, ])
        n_sub_unique <- nrow(sub_unique)
        nA_sub       <- length(unique(sub_unique$X_A))
        nB_sub       <- length(unique(sub_unique$X_B))
        max_poss     <- nA_sub * nB_sub
        dens_sub     <- n_sub_unique / max_poss
        
        sim_res[[j]] <- data.frame(
          Simulation  = i,
          Size        = n,
          UniqueEdges = n_sub_unique,
          nA          = nA_sub,
          nB          = nB_sub,
          MaxPossible = max_poss,
          Density     = dens_sub
        )
      }
      
      message(sprintf("[sim %d/%d] DONE", i, nsim))
      do.call(rbind, sim_res)
    },
    mc.cores      = ncores,
    mc.preschedule= FALSE  # dispatch tasks as workers free up
  )
  
  message("Combining all simulation results...")
  result_df <- do.call(rbind, results_list)
  message("All done.")
  
  return(result_df)
}

dens_PY25 <- dens(sim_edgesPY552525)
dens_PY5 <- dens(sim_edgesPY5555)
dens_PY75 <- dens(sim_edgesPY557575)

save(dens_PY25, file = here(dens_folder, "dens_PY_0.25.Rdata"))
save(dens_PY5, file = here(dens_folder, "dens_PY_0.5.Rdata"))
save(dens_PY75, file = here(dens_folder, "dens_PY_0.75.Rdata"))