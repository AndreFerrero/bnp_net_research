# 0. Load dependencies
if (!requireNamespace("pbmcapply", quietly = TRUE)) {
  install.packages("pbmcapply")
}
library(pbmcapply)
source("code/setup.R")

# 1. Master settings
set.seed(123)
nsim         <- 200
N            <- 1e9
alpha_values <- c(0.5, 1, 1.5, 2, 3, 5, 20)
out_dir      <- "sims"
numCores     <- 80

# 2. Grand‐total timer
total_start <- Sys.time()
message("=== ALL SIMULATIONS STARTING at ", total_start, " ===\n")

# 3. Loop over α’s
for (i in seq_along(alpha_values)) {
  alpha_val <- alpha_values[i]
  a_str     <- gsub("\\.", "", as.character(alpha_val))
  suffix    <- paste0(a_str, a_str)
  sim_name  <- paste0("sim_edgesDP", suffix)
  file_out  <- file.path(out_dir, paste0("dp", suffix, ".Rdata"))

  # log start
  message(sprintf("  [%d/%d] START α = %-4s  at %s",
                  i, length(alpha_values),
                  alpha_val, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()

  # run nsim simulations in parallel with progress bar
  sims_info <- pbmcapply::pbmclapply(
    X         = seq_len(nsim),
    mc.cores  = numCores,
    mc.preschedule = FALSE,
    FUN       = function(j) {
      t_j0 <- Sys.time()
      net  <- sample_net(N,
                         alpha = c(alpha_val, alpha_val),
                         sigma = c(0, 0))
      elapsed <- as.numeric(difftime(Sys.time(), t_j0, units = "secs"))
      n_edges <- if (is.data.frame(net)) nrow(net) else NA_integer_
      list(net = net, elapsed = elapsed, n_edges = n_edges)
    }
  )

  # Extract nets + print sim summaries
  nets <- vector("list", nsim)
  for (j in seq_len(nsim)) {
    info      <- sims_info[[j]]
    nets[[j]] <- info$net
    cat(sprintf("    [α=%s] sim %3d/%3d DONE in %5.2f secs → %d edges\n",
                alpha_val, j, nsim,
                info$elapsed,
                info$n_edges))
  }

  # Save result
  assign(sim_name, nets, envir = .GlobalEnv)
  save(list = sim_name, file = file_out)

  # log batch end
  elapsed_batch <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  [%d/%d] DONE  α = %-4s  elapsed = %5.1f secs → %s\n",
                  i, length(alpha_values),
                  alpha_val, elapsed_batch,
                  basename(file_out)))
}

# 4. Done
total_elapsed <- as.numeric(difftime(Sys.time(), total_start, units = "mins"))
message("=== ALL SIMULATIONS FINISHED at ", Sys.time(),
        " (total = ", round(total_elapsed, 1), " mins) ===")

