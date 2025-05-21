# Full R script: one‐time cluster + live progress bars per α

# 0. Install pbapply if you haven’t already
# install.packages("pbapply")

# 1. Load required packages
library(parallel)
library(pbapply)

# 2. Master settings
set.seed(123)                              # reproducibility
nsim         <- 50                       # sims per α
N            <- 1e6                        # edges per sim
alpha_values <- c(0.5, 1, 2, 3 ,5, 20)     # α’s to loop over
out_dir      <- "C:/Users/anferrer/OneDrive - Alma Mater Studiorum Università di Bologna/University/Internship/research/sims"

# 3. Create cluster once & export static objects
numCores <- 6
cl       <- makeCluster(numCores)
clusterExport(cl, varlist = c("sample_net", "N", "py_predictive_draw3"),
              envir = environment())

# 4. Grand‐total timer
total_start <- Sys.time()
message("=== ALL SIMULATIONS STARTING at ", total_start, " ===\n")

# 5. Loop over α’s
for (i in seq_along(alpha_values)) {
  alpha_val <- alpha_values[i]
  a_str     <- gsub("\\.", "", as.character(alpha_val))
  suffix    <- paste0(a_str, a_str)               # e.g. 05 → "0505"
  sim_name  <- paste0("sim_edgesDP", suffix)      # in‐workspace name
  file_out  <- file.path(out_dir, paste0("dp", suffix, ".Rdata"))
  
  # high‐level start
  message(sprintf("  [%d/%d] START α = %-4s  at %s",
                  i, length(alpha_values),
                  alpha_val, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()
  
  # tell workers the new alpha
  clusterExport(cl, varlist = "alpha_val", envir = environment())
  
  # set up a console timer/progress bar for the next pblapply
  pboptions(type = "timer", style = 3)
  
  # run sims with live progress bar
  sims_info <- pblapply(
    X   = seq_len(nsim),
    cl  = cl,
    FUN = function(j) {
      t_j0    <- Sys.time()
      net     <- sample_net(N,
                            alpha = c(alpha_val, alpha_val),
                            sigma = c(0, 0))
      elapsed <- as.numeric(difftime(Sys.time(), t_j0, units = "secs"))
      n_edges <- if (is.data.frame(net)) nrow(net) else NA_integer_
      list(net     = net,
           elapsed = elapsed,
           n_edges = n_edges)
    }
  )
  
  # extract nets & print a summary line per sim
  nets <- vector("list", nsim)
  for (j in seq_len(nsim)) {
    info      <- sims_info[[j]]
    nets[[j]] <- info$net
    cat(sprintf("    [α=%s] sim %3d/%3d DONE in %5.2f secs → %d edges\n",
                alpha_val, j, nsim,
                info$elapsed,
                info$n_edges))
  }
  
  # assign + save
  assign(sim_name, nets, envir = .GlobalEnv)
  save(list = sim_name, file = file_out)
  
  # high‐level finish for this α
  elapsed_batch <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  [%d/%d] DONE  α = %-4s  elapsed = %5.1f secs → %s\n",
                  i, length(alpha_values),
                  alpha_val, elapsed_batch,
                  basename(file_out)))
}

# 6. Tear down & grand‐total finish
stopCluster(cl)
total_elapsed <- as.numeric(difftime(Sys.time(), total_start, units = "mins"))
message("=== ALL SIMULATIONS FINISHED at ", Sys.time(),
        " (total = ", round(total_elapsed,1), " mins) ===")
