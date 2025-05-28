# Load necessary libraries
library(foreach)
library(doParallel)
library(dplyr)
library(utils)
library(here)

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
  xA <- py_sample(N, alpha[1], sigma[1])
  xB <- py_sample(N, alpha[2], sigma[2])
  edges <- data.frame(X_A = xA$x, X_B = xB$x)
  list(edges = edges, xA = xA, xB = xB)
}

compute_density <- function(net) {
  unique_edges <- nrow(unique(net$edges))
  unique_edges / (net$xA$K * net$xB$K)
}

# SIMULATION PARAMETERS --------------------------------------------------
sizes <- 2^(4:17)
n_rep <- 50
alpha <- c(5, 5)
sigma_list <- list(c(0.25, 0.25), c(0.5, 0.5), c(0.75, 0.75))

# Set root results directory using here
root_dir <- here("py_dens", "sim_results")
if (!dir.exists(root_dir)) dir.create(root_dir, recursive = TRUE)

# Parallel backend setup
n_cores <- min(n_rep, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# MAIN SIMULATION ---------------------------------------------------------
all_summaries <- list()
total_steps <- length(sigma_list) * length(sizes)
step_counter <- 0
cat("Starting simulation with", total_steps, "steps...\n")
pb <- txtProgressBar(min = 0, max = total_steps, style = 3)

for (s in sigma_list) {
  sigmaA <- s[1]
  sigmaB <- s[2]
  sigma_str <- paste0(sigmaA, "_", sigmaB)

  for (N in sizes) {
    # prepare directories under root_dir
    base_dir <- file.path(root_dir, paste0("sigma_", sigma_str), paste0("size_", N))
    dir.create(file.path(base_dir, "networks"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(base_dir, "densities"), recursive = TRUE, showWarnings = FALSE)

    # parallel over replicates
    results_N <- foreach(
      rep = seq_len(n_rep), .combine = rbind,
      .export = c("py_sample", "sample_net", "compute_density")
    ) %dopar% {
      net <- sample_net(N, alpha, c(sigmaA, sigmaB))

      net_file <- file.path(base_dir, "networks", paste0("network_rep_", rep, ".rds"))
      saveRDS(net, net_file)

      dens <- compute_density(net)
      dens_file <- file.path(base_dir, "densities", paste0("density_rep_", rep, ".rds"))
      saveRDS(dens, dens_file)

      data.frame(
        sigmaA  = sigmaA,
        sigmaB  = sigmaB,
        size    = N,
        rep     = rep,
        density = dens
      )
    }

    all_summaries[[paste0(sigma_str, "_", N)]] <- results_N
    step_counter <- step_counter + 1
    setTxtProgressBar(pb, step_counter)
  }
}
close(pb)
stopCluster(cl)

# AGGREGATE STATISTICS ----------------------------------------------------
results_df <- bind_rows(all_summaries)
summary_stats <- results_df %>%
  group_by(sigmaA, sigmaB, size) %>%
  summarise(
    mean_density = mean(density),
    sd_density   = sd(density),
    se           = sd_density / sqrt(n()),
    ci_lower     = mean_density - qt(0.975, n() - 1) * se,
    ci_upper     = mean_density + qt(0.975, n() - 1) * se,
    .groups      = "drop"
  )

# Save summary under root_dir
saveRDS(summary_stats, file.path(root_dir, "density_summary.rds"))
write.csv(summary_stats, file.path(root_dir, "density_summary.csv"), row.names = FALSE)

# PLOTTING ---------------------------------------------------------------
summary_stats <- read.csv(here("res", "density_PY", "new_sim", "density_summary.csv"))
library(ggplot2)

library(ggplot2)

# Compute log2(size) as a separate variable for clarity
summary_stats$log2_size <- log2(summary_stats$size)

p_dens <- ggplot(summary_stats, aes(x = log2_size, y = mean_density, color = factor(sigmaA))) +
  geom_line(size = 1) +
  geom_point(aes(shape = factor(sigmaA)), size = 3, stroke = 0.8) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(sigmaA)), alpha = 0.2, color = NA) +
  scale_shape_manual(values = rep(15, length(unique(summary_stats$sigmaA)))) +
  scale_x_continuous(
    breaks = log2(unique(summary_stats$size))
  ) +
  labs(
    title = "Bipartite Network Density vs. Size",
    x = expression(log[2](n)),
    y = "Density",
    color = expression(sigma),
    fill = expression(sigma),
    shape = expression(sigma)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


ggsave(p_dens, filename = here("res", "pics", "density_analysis", "py", "log2_dens_summary.pdf"))
# End of script
