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
  unique_edges / (as.numeric(net$xA$K) * as.numeric(net$xB$K))
}

# SIMULATION PARAMETERS --------------------------------------------------
sizes <- 2^(4:17)
n_rep <- 50
alpha <- c(5, 5)
sigma_vals <- seq(0.1, 0.9, 0.1)
sigma_list <- lapply(sigma_vals, function(x) c(x, x))

# Set root results directory using here
root_dir <- here("py_dens", "rate_results")
if (!dir.exists(root_dir)) dir.create(root_dir, recursive = TRUE)

# Parallel backend setup
n_cores <- min(n_rep, parallel::detectCores())
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
    # parallel over replicates
    results_N <- foreach(
      rep = seq_len(n_rep), .combine = rbind,
      .export = c("py_sample", "sample_net", "compute_density")
    ) %dopar% {
      net <- sample_net(N, alpha, c(sigmaA, sigmaB))
      dens <- compute_density(net)

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
dens_summary <- results_df |>
  group_by(sigmaA, sigmaB, size) |>
  summarise(
    mean_density = mean(density),
    sd_density   = sd(density),
    se           = sd_density / sqrt(n()),
    ci_lower     = mean_density - qt(0.975, n() - 1) * se,
    ci_upper     = mean_density + qt(0.975, n() - 1) * se,
    .groups      = "drop"
  )

# Add log2-transformed variables
dens_summary <- dens_summary |>
  mutate(
    log2_size = log2(size),
    log2_density = log2(mean_density),
    log2_ci_lower = log2(ci_lower),
    log2_ci_upper = log2(ci_upper),
    sigma = paste0("sigma = ", sigmaA) # For grouping/plotting
  )

slopes_df <- dens_summary |>
  group_by(sigmaA, sigmaB) |>
  do({
    model <- lm(log2_density ~ log2_size, data = .)
    slope <- coef(model)[["log2_size"]]
    intercept <- coef(model)[["(Intercept)"]]
    data.frame(slope = slope, intercept = intercept)
  }) |>
  ungroup()

slopes_df <- slopes_df |>
  mutate(sigma = sigmaA)

write.csv(slopes_df, file = file.path(root_dir, "slopes_summary.csv"), row.names = FALSE)
write.csv(dens_summary, file = file.path(root_dir, "dens_summary.csv"), row.names = FALSE)


slopes_df <- read.csv(here("res", "density_PY", "new_sim", "rate_results", "slopes_summary.csv"))
dens_summary <- read.csv(here("res", "density_PY", "new_sim", "rate_results", "dens_summary.csv"))

library(ggplot2)

ggplot(slopes_df, aes(x = sigma, y = slope)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  scale_x_continuous(
    breaks = seq(0.1, 0.9, 0.1)
  ) +
  theme_minimal(base_size = 14)

ggplot(dens_summary, aes(x = log2_size, y = mean_density, color = factor(sigmaA))) +
  geom_line(size = 0.1) +
  geom_point(aes(shape = factor(sigmaA)), size = 0.5, stroke = 0.8) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = factor(sigmaA)), alpha = 0.2, color = NA) +
  scale_x_continuous(
    breaks = log2(c(100, 1000, 10000, 100000)),
    labels = c(
      expression(10^2),
      expression(10^3),
      expression(10^4),
      expression(10^5)
    )
  ) +
  labs(
    x = "n",
    y = "density",
    color = expression(sigma),
    fill = expression(sigma),
    shape = expression(sigma)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggplot(dens_summary, aes(x = log2_size, y = log2_density, color = factor(sigmaA))) +
  geom_line(size = 1) +
  geom_point(aes(shape = factor(sigmaA)), size = 3, stroke = 0.8) +
  geom_ribbon(aes(ymin = log2_ci_lower, ymax = log2_ci_upper, fill = factor(sigmaA)), alpha = 0.2, color = NA) +

  # X-axis: log2 of size with 10^j labels
  scale_x_continuous(
    breaks = log2(c(100, 1000, 10000, 100000)),
    labels = c(expression(10^2), expression(10^3), expression(10^4), expression(10^5))
  ) +

  # Y-axis: log2 of density with parsed labels if you want similar scientific style
  scale_y_continuous(
    breaks = pretty(dens_summary$log2_density),
    labels = function(x) parse(text = paste0("2^", x))
  ) +
  labs(
    x = "n", # Keep the label simple
    y = "density", # Since log2 scale is encoded in ticks
    color = expression(sigma),
    fill = expression(sigma),
    shape = expression(sigma)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
