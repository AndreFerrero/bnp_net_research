library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(tidyverse)
library(posterior)

# Directories
poll_dir <- here("pollinators/mpl006")
w_repl_dir <- here(poll_dir, "w_repl")
plots_dir <- here(w_repl_dir, "beta_01_1_plots")
fits_dir <- here(w_repl_dir, "fits")

dir.create(plots_dir, showWarnings = FALSE)

# Read original full‐edge list once
raw_mat <- read.csv(here(poll_dir, "M_PL_006.csv"),
  check.names = FALSE, row.names = 1
)
full_edges <- raw_mat %>%
  rownames_to_column("plant") %>%
  pivot_longer(-plant, names_to = "pollinator", values_to = "weight") %>%
  filter(weight > 0)

total_weight <- sum(full_edges$weight)

# Your sample sizes (not percentages now but exact counts)
sample_sizes <- c(10000, 100000, 1000000)
delta <- 0.1

# Containers to collect posterior draws and Bayes factors
sigma_summaries <- tibble()
bf_summary <- tibble()

for (n in sample_sizes) {
  cat("Processing sample size", n, "...\n")

  # Adjust file name accordingly
  file_name <- paste0("beta_01_1_fit_", sprintf("%d", n), ".Rdata")
  load(here(fits_dir, file_name))

  # diagnostics: trace, acf, dens overlay for parameters
  pars <- c("alpha_A", "alpha_B", "sigma_A", "sigma_B")

  trace_plot <- mcmc_trace(fit, pars = pars) +
    ggtitle(paste0(n, " samples trace")) +
    theme(legend.position = "top")
  ggsave(here(plots_dir, paste0("trace_", n, ".pdf")), trace_plot)

  acf_plot <- mcmc_acf_bar(fit, pars = pars) +
    ggtitle(paste0(n, " samples ACF")) +
    theme(legend.position = "top")
  ggsave(here(plots_dir, paste0("acf_", n, ".pdf")), acf_plot)

  dens_plot <- mcmc_dens_overlay(fit, pars = pars) +
    ggtitle(paste0(n, " samples density overlay")) +
    theme(legend.position = "top")
  ggsave(here(plots_dir, paste0("dens_overlay_", n, ".pdf")), dens_plot)

  # Extract draws (no PPC)
  draws_df <- as_draws_df(fit)

  # Bayes factor for sigma_A < delta
  post0 <- mean(draws_df$sigma_A < delta)
  post1 <- mean(draws_df$sigma_A > delta)
  prior0 <- pbeta(delta, 0.05, 1)
  prior1 <- 1 - prior0
  bf <- (post0 / post1) * (prior1 / prior0)
  log10_bf <- log10(bf)

  # Save posterior sigma draws and BF summary
  sigma_summaries <- bind_rows(
    sigma_summaries,
    tibble(
      sample_size = n,
      sigma_A = draws_df$sigma_A,
      sigma_B = draws_df$sigma_B
    )
  )
  bf_summary <- bind_rows(
    bf_summary,
    tibble(
      sample_size = n,
      BF = bf,
      log10BF = log10_bf
    )
  )
}

# ---- Plot posterior densities of sigma_A and sigma_B across sample sizes ----
sigma_long <- sigma_summaries %>%
  pivot_longer(c(sigma_A, sigma_B), names_to = "parameter", values_to = "value")

dens_all <- ggplot(sigma_long, aes(x = value, color = factor(sample_size))) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  labs(
    color = "Sample size",
    title = "Posterior densities of sigma_A and sigma_B across sample sizes"
  ) +
  theme_minimal()
ggsave(here(plots_dir, "sigma_densities_all.pdf"), dens_all, width = 8, height = 4)

# ---- Plot Bayes factor vs. sample size ----
bf_plot <- ggplot(bf_summary, aes(x = sample_size, y = log10BF)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(
    x = "Sample size",
    y = "log10 Bayes Factor",
    title = bquote("Bayes factor for " ~ sigma[A] < .(delta))
  ) +
  theme_minimal()
ggsave(here(plots_dir, "bayes_factor_vs_sample_size.pdf"), bf_plot, width = 6, height = 4)
