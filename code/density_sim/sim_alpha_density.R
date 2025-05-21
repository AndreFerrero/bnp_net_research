# density_simulation_with_rcpp.R

# --- 0. Load dependencies and compile C++ implementation of sample_net ---
if (!requireNamespace("Rcpp", quietly = TRUE))   install.packages("Rcpp")
if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply")
if (!requireNamespace("pbapply", quietly = TRUE))    install.packages("pbapply")

library(pbmcapply)
library(pbapply)

source("code/1_py_predictive_draw3.R")
source("code/3_sample_net.R")

# --- 1. Master settings ---
set.seed(123)
nsim         <- 500
N            <- 1e7
alpha_values <- c(0.5, 1, 1.5, 2, 3, 5, 20)
thresholds   <- c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7)
out_dir      <- "sims500"
numCores     <- parallel::detectCores()

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 2. Grand-total timer ---
total_start <- Sys.time()
message("=== ALL SIMULATIONS STARTING at ", total_start, " ===\n")

# --- 3. Simulate networks for each alpha ---
sims_by_alpha <- list()

for (i in seq_along(alpha_values)) {
  alpha_val <- alpha_values[i]
  a_str     <- gsub("\\.", "", as.character(alpha_val))
  suffix    <- paste0(a_str, a_str)
  sim_name  <- paste0("sim_edgesDP", suffix)
  file_out  <- file.path(out_dir, paste0("dp", suffix, ".Rdata"))
  
  message(sprintf("  [%d/%d] START α = %-4s  at %s",
                  i, length(alpha_values),
                  alpha_val, format(Sys.time(), "%H:%M:%S")))
  t0 <- Sys.time()
  
  sims_info <- pbmclapply(
    seq_len(nsim), mc.cores = numCores, mc.preschedule = FALSE,
    FUN = function(j) {
      t_j0 <- Sys.time()
      net  <- sample_net(N, c(alpha_val, alpha_val), c(0, 0))
      elapsed <- as.numeric(difftime(Sys.time(), t_j0, units = "secs"))
      n_edges <- if (is.data.frame(net)) nrow(net) else NA_integer_
      list(net = net, elapsed = elapsed, n_edges = n_edges)
    }
  )
  
  nets <- vector("list", nsim)
  for (j in seq_len(nsim)) {
    info      <- sims_info[[j]]
    nets[[j]] <- info$net
    cat(sprintf("    [α=%s] sim %3d/%3d DONE in %5.2f secs → %d edges\n",
                alpha_val, j, nsim,
                info$elapsed,
                info$n_edges))
  }
  
  assign(sim_name, nets, envir = .GlobalEnv)
  save(list = sim_name, file = file_out)
  sims_by_alpha[[as.character(alpha_val)]] <- nets
  
  elapsed_batch <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  [%d/%d] DONE  α = %-4s  elapsed = %5.1f secs → %s\n",
                  i, length(alpha_values),
                  alpha_val, elapsed_batch,
                  basename(file_out)))
}

total_elapsed <- as.numeric(difftime(Sys.time(), total_start, units = "mins"))
message("=== ALL SIMULATIONS FINISHED at ", Sys.time(),
        " (total = ", round(total_elapsed, 1), " mins) ===")

# --- 4. Density analysis with progress messages ---
avg_density_by_alpha <- function(sims_by_alpha, thresholds, ncores) {
  jobs <- expand.grid(alpha = names(sims_by_alpha),
                      threshold = thresholds,
                      stringsAsFactors = FALSE)
  total_jobs <- nrow(jobs)
  message(sprintf("=== Starting density computation: %d tasks ===", total_jobs))
  
  results_list <- pbmclapply(
    seq_len(total_jobs),
    mc.cores       = min(total_jobs, ncores),
    mc.preschedule = FALSE,
    FUN = function(i) {
      a   <- jobs$alpha[i]
      th  <- jobs$threshold[i]
      # progress message
      message(sprintf("[Density %d/%d] α=%s, threshold=%s", 
                      i, total_jobs, a, format(th, scientific=FALSE)))
      sims <- sims_by_alpha[[a]]
      dens <- vapply(sims, function(net) {
        sb <- unique(net[1:th, , drop = FALSE])
        ne <- nrow(sb)
        nA <- length(unique(sb$X_A))
        nB <- length(unique(sb$X_B))
        ne / (nA * nB)
      }, numeric(1))
      data.frame(alpha = as.numeric(a),
                 threshold = th,
                 mean_density = mean(dens))
    }
  )
  do.call(rbind, results_list)
}

# --- 5. Run density analysis and save ---
alpha_dens <- avg_density_by_alpha(sims_by_alpha, thresholds,
                                   ncores = numCores)
save(alpha_dens, file = file.path(out_dir, "50alpha_density_results.Rdata"))
message("=== Density results saved to ", file.path(out_dir), " ===")
