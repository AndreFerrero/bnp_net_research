# Load necessary libraries
library(foreach)
library(doParallel)
library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")

# FUNCTIONS ---------------------------------------------------------------
py_sample <- function(N, alpha = 5, sigma = 0) {
  x <- numeric(N)
  x[1] <- 1
  counts <- numeric(N)
  counts[1] <- 1
  active_nodes <- c(1)

  for (n in 2:N) {
    k <- length(active_nodes)
    prob <- c(counts[active_nodes] - sigma, alpha + sigma * k)
    x[n] <- sample(c(active_nodes, k + 1), 1, prob = prob)
    counts[x[n]] <- counts[x[n]] + 1
    if (counts[x[n]] == 1) active_nodes <- c(active_nodes, x[n])
  }

  list(
    x = x,
    counts = counts,
    active_nodes = active_nodes,
    active_counts = counts[active_nodes],
    K = length(active_nodes)
  )
}

sample_net <- function(N, alpha, sigma = c(0, 0)) {
  alphaA <- alpha[1]
  alphaB <- alpha[2]
  sigmaA <- sigma[1]
  sigmaB <- sigma[2]
  xA <- py_sample(N, alphaA, sigmaA)
  xB <- py_sample(N, alphaB, sigmaB)
  edges <- data.frame(X_A = xA$x, X_B = xB$x)
  list(edges = edges, xA = xA, xB = xB)
}

compute_density <- function(net) {
  unique_edges <- nrow(unique(net$edges))
  K_A <- net$xA$K
  K_B <- net$xB$K
  unique_edges / (K_A * K_B)
}

# SIMULATION PARAMETERS --------------------------------------------------
sizes <- c(1e3, 1e4, 1e5, 1e6)
n_rep <- 50
alpha <- c(5, 5)
sigma_list <- list(c(0.25, 0.25), c(0.5, 0.5), c(0.75, 0.75))

# Parallel backend setup
n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

# MAIN SIMULATION ---------------------------------------------------------
all_summaries <- list()

progressr::with_progress({
  p <- progressor(steps = length(sigma_list) * length(sizes))

  for (s in sigma_list) {
    sigmaA <- s[1]
    sigmaB <- s[2]
    sigma_str <- paste0(sigmaA, "_", sigmaB)

    for (N in sizes) {
      # prepare directories
      base_dir <- file.path("sim_results", paste0("sigma_", sigma_str), paste0("size_", N))
      dir.create(file.path(base_dir, "networks"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(base_dir, "densities"), recursive = TRUE, showWarnings = FALSE)

      # parallel over replicates for this size
      results_N <- foreach(
        rep = seq_len(n_rep), .combine = rbind,
        .packages = c("progressr"),
        .export = c("sample_net", "compute_density")
      ) %dopar% {
        progressr::handlers("txtprogressbar")
        progressr::with_progress({
          p_inner <- progressor(steps = 1)

          # simulate
          net <- sample_net(N, alpha, c(sigmaA, sigmaB))

          # save network
          net_file <- file.path(base_dir, "networks", paste0("network_rep_", rep, ".rds"))
          saveRDS(net, net_file)

          # compute & save density
          dens <- compute_density(net)
          dens_file <- file.path(base_dir, "densities", paste0("density_rep_", rep, ".rds"))
          saveRDS(dens, dens_file)

          # mark progress
          p_inner()

          # return record
          data.frame(
            sigmaA = sigmaA,
            sigmaB = sigmaB,
            size = N,
            rep = rep,
            density = dens
          )
        })
      }

      all_summaries[[paste0(sigma_str, "_", N)]] <- results_N
      p(sprintf("Completed sigma=(%.2f, %.2f), N=%d", sigmaA, sigmaB, N))
    }
  }
})

# stop cluster
stopCluster(cl)

# AGGREGATE STATISTICS ----------------------------------------------------
library(dplyr)
results_df <- bind_rows(all_summaries)
summary_stats <- results_df %>%
  group_by(sigmaA, sigmaB, size) %>%
  summarise(
    mean_density = mean(density),
    sd_density = sd(density),
    se = sd_density / sqrt(n()),
    ci_lower = mean_density - qt(0.975, n() - 1) * se,
    ci_upper = mean_density + qt(0.975, n() - 1) * se,
    .groups = "drop"
  )

# save summary
dir.create("sim_results", showWarnings = FALSE)
saveRDS(summary_stats, file.path("sim_results", "density_summary.rds"))
write.csv(summary_stats, file.path("sim_results", "density_summary.csv"), row.names = FALSE)
# End of script
