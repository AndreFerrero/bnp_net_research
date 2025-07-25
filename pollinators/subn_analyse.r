library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(tidyverse)
library(posterior)

# Set up directories
poll_dir <- here("pollinators/mpl025")
wo_repl_dir <- here(poll_dir, "wo_repl")
plots_dir <- here(wo_repl_dir, "unif_plots")
fits_dir <- here(wo_repl_dir, "unif_fit")

dir.create(plots_dir, showWarnings = FALSE)

# Load full data
raw_mat <- read.csv(here(poll_dir, "M_PL_025.csv"),
  check.names = FALSE, row.names = 1
)

full_edges <- raw_mat %>%
  rownames_to_column("plant") %>%
  pivot_longer(-plant, names_to = "pollinator", values_to = "weight") %>%
  filter(weight > 0)

total_weight <- sum(full_edges$weight)

# Compute observed density
compute_d_obs <- function(edges, subn) {
  tickets <- rep(seq_len(nrow(edges)), times = edges$weight)
  sel <- sample(tickets, size = subn, replace = FALSE)
  df <- as.data.frame(table(sel), stringsAsFactors = FALSE) %>%
    transmute(row_index = as.integer(sel), weight = as.integer(Freq))
  sub_e <- edges %>%
    slice(df$row_index) %>%
    mutate(weight = df$weight)
  sub_pl <- sub_e %>%
    group_by(plant) %>%
    summarise(counts = sum(weight))
  sub_po <- sub_e %>%
    group_by(pollinator) %>%
    summarise(counts = sum(weight))
  list(
    d_obs = nrow(sub_e) / (nrow(sub_pl) * nrow(sub_po)),
    e_obs = sum(sub_e$weight)
  )
}

# Parameters to iterate
percentages <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
delta_values <- c(0.1, 0.05, 0.01)

# Storage
sigma_summaries <- tibble()
bf_all_deltas <- tibble()

# Main loop over fits
for (p in percentages) {
  pct_label <- sprintf("%03d", round(100 * p))
  cat("Processing", pct_label, "…\n")
  load(here(fits_dir, paste0("unif_ppc_fit_", pct_label, "pct.Rdata")))

  # Recompute d_obs
  ss <- compute_d_obs(full_edges, round(total_weight * p))

  # Diagnostics
  pars <- c("alpha_A", "alpha_B", "sigma_A", "sigma_B")

  trace_plot <- mcmc_trace(fit, pars = pars) +
    ggtitle(paste0(pct_label, "% trace")) +
    theme(legend.position = "top")
  ggsave(here(plots_dir, paste0("trace_", pct_label, "pct.pdf")), trace_plot)

  acf_plot <- mcmc_acf_bar(fit, pars = pars) +
    ggtitle(paste0(pct_label, "% ACF")) +
    theme(legend.position = "top")
  ggsave(here(plots_dir, paste0("acf_", pct_label, "pct.pdf")), acf_plot)

  dens_plot <- mcmc_dens_overlay(fit, pars = pars) +
    ggtitle(paste0(pct_label, "% dens overlay")) +
    theme(legend.position = "top")
  ggsave(here(plots_dir, paste0("dens_overlay_", pct_label, "pct.pdf")), dens_plot)

  # Posterior predictive check
  draws_df <- as_draws_df(fit)
  d_ppc <- draws_df$density_ppc
  ppc_hist <- ggplot(data.frame(density = d_ppc), aes(x = density)) +
    geom_histogram(bins = 30, color = "black", fill = "lightblue") +
    geom_vline(xintercept = ss$d_obs, color = "red", size = 1) +
    labs(
      title = paste0(pct_label, "% PPC density"),
      x = expression(d == e / (K[A] * K[B]))
    ) +
    theme_minimal()
  ggsave(here(plots_dir, paste0("ppc_density_", pct_label, "pct.pdf")), ppc_hist)

  # Posterior storage
  sigma_summaries <- bind_rows(
    sigma_summaries,
    tibble(
      pct        = 100 * p,
      sigma_A    = draws_df$sigma_A,
      sigma_B    = draws_df$sigma_B
    )
  )

  # Compute Bayes Factors for all deltas
  for (delta in delta_values) {
    # Posterior probabilities
    post0A <- mean(draws_df$sigma_A < delta)
    post1A <- mean(draws_df$sigma_A > delta)
    post0B <- mean(draws_df$sigma_B < delta)
    post1B <- mean(draws_df$sigma_B > delta)

    # Prior (same for both groups)
    prior0 <- pbeta(delta, 1, 1)
    prior1 <- 1 - prior0

    # Bayes Factors
    bf_A <- (post0A / post1A) * (prior1 / prior0)
    log10_bf_A <- log10(bf_A)

    bf_B <- (post0B / post1B) * (prior1 / prior0)
    log10_bf_B <- log10(bf_B)

    # Store all
    bf_all_deltas <- bind_rows(
      bf_all_deltas,
      tibble(
        pct        = 100 * p,
        delta      = delta,
        BF_A       = bf_A,
        log10BF_A  = log10_bf_A,
        BF_B       = bf_B,
        log10BF_B  = log10_bf_B
      )
    )
  }
}

# Plot posterior densities
library(viridis)

sigma_long <- sigma_summaries %>%
  pivot_longer(c(sigma_A, sigma_B), names_to = "parameter", values_to = "value")

parameter_labels <- c(
  sigma_A = "sigma[A]",
  sigma_B = "sigma[B]"
)

dens_all <- ggplot(sigma_long, aes(x = value, color = factor(pct))) +
  geom_density(key_glyph = "path") +
  facet_wrap(
    ~parameter,
    scales = "free",
    labeller = labeller(parameter = as_labeller(parameter_labels, label_parsed))
  ) +
  scale_color_viridis_d(
    option = "B",
    end = 0.9,
    direction = -1,
    guide = guide_legend(nrow = 1, override.aes = list(linetype = 1, shape = NA, fill = NA))
  ) +
  labs(
    color = "sample size (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.y = element_blank(),
    axis.title.x = element_blank()
  )


ggsave(here(plots_dir, "sigma_densities_all_viridisB.pdf"), dens_all, width = 8, height = 4)

# ---- Separate Bayes Factor plots with fixed y-axis ticks and limits ----
y_ticks <- c(0, 0.5, 1, 2)

for (delta in delta_values) {
  df_delta <- bf_all_deltas %>% filter(delta == !!delta)
  delta_label <- str_replace(format(delta, nsmall = 2), "\\.", "")

  bf_plot <- ggplot(df_delta, aes(x = pct, y = log10BF)) +
    geom_line() +
    geom_point() +
    labs(
      x = "sample size (%)",
      y = "log10 Bayes Factor",
      title = bquote("Bayes factor for " ~ sigma[A] < .(delta))
    ) +
    scale_x_continuous(
      breaks = unique(df_delta$pct),
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      breaks = y_ticks,
      limits = c(0, 2), # force showing tick at 2
      expand = c(0.01, 0)
    ) +
    theme_minimal()

  ggsave(
    here(plots_dir, paste0("bayes_factor_delta_", delta_label, ".pdf")),
    bf_plot,
    width = 6, height = 4
  )
}

# ---- Combined Bayes Factor plot with fixed y-axis limits ----
library(viridis)
y_ticks_A <- c(0, 0.5, 1, 2)

# Combined Bayes Factor plot — Group A
bf_combined_plot_A <- ggplot(bf_all_deltas, aes(x = pct, y = log10BF_A)) +
  geom_line(aes(linetype = factor(delta), group = factor(delta)), color = "black") +
  geom_point(aes(color = factor(pct))) +
  scale_color_viridis_d(option = "B", end = 0.9, direction = -1, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = expression(delta)) +
  labs(
    x = "sample size (%)",
    y = expression(log[10](BF)),
    title = expression("Bayes factor for" ~ sigma[A])
  ) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = Inf, alpha = 0.1, fill = "gray50") +
  scale_x_continuous(breaks = unique(bf_all_deltas$pct), expand = c(0.01, 0)) +
  scale_y_continuous(breaks = y_ticks_A, limits = c(0, 2.5), expand = c(0.01, 0)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank()
  )

ggsave(
  here(plots_dir, "bayes_factor_all_deltas_groupA.pdf"),
  bf_combined_plot_A,
  width = 7, height = 4
)

y_ticks_B <- c(-0.5, 0, 0.5, 1)

# Combined Bayes Factor plot — Group B
bf_combined_plot_B <- ggplot(bf_all_deltas, aes(x = pct, y = log10BF_B)) +
  geom_line(aes(linetype = factor(delta), group = factor(delta)), color = "black") +
  geom_point(aes(color = factor(pct))) +
  scale_color_viridis_d(option = "B", end = 0.9, direction = -1, guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = expression(delta)) +
  labs(
    x = "sample size (%)",
    y = expression(log[10](BF)),
    title = expression("Bayes factor for" ~ sigma[B])
  ) +
  scale_x_continuous(breaks = unique(bf_all_deltas$pct), expand = c(0.01, 0)) +
  scale_y_continuous(breaks = y_ticks_B, limits = c(-0.5, 1), expand = c(0.01, 0)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank()
  )

ggsave(
  here(plots_dir, "bayes_factor_all_deltas_groupB.pdf"),
  bf_combined_plot_B,
  width = 7, height = 4
)
