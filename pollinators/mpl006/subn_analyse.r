library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(tidyverse)
library(posterior)

# where things live
poll_dir <- here("pollinators/mpl006")
wo_repl_dir <- here(poll_dir, "wo_repl")
plots_dir <- here(wo_repl_dir, "hyper_beta_plots")
fits_dir <- here(wo_repl_dir, "fits")

dir.create(plots_dir)

# read original full‐edge list once
raw_mat <- read.csv(here(poll_dir, "M_PL_006.csv"),
  check.names = FALSE, row.names = 1
)
full_edges <- raw_mat %>%
  rownames_to_column("plant") %>%
  pivot_longer(-plant, names_to = "pollinator", values_to = "weight") %>%
  filter(weight > 0)

total_weight <- sum(full_edges$weight)

# helper to recompute d_obs at a given subn
compute_d_obs <- function(edges, subn) {
  # identical to your sub() minus summarise
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

# percentages to loop
percentages <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)
delta <- 0.1

# containers for summary across fits
sigma_summaries <- tibble()
bf_summary <- tibble()

for (p in percentages) {
  pct_label <- sprintf("%03d", round(100 * p))
  cat("Processing", pct_label, "…\n")
  load(here(fits_dir, paste0("hyper_ppc_fit_", pct_label, "pct.Rdata")))

  # recompute d_obs & e_obs
  ss <- compute_d_obs(full_edges, round(total_weight * p))

  # diagnostics: trace, acf, dens overlay for 4 parameters
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

  # PPC on network density
  # draws_df <- as_draws_df(fit)
  # d_ppc <- draws_df$density_ppc
  # ppc_hist <- ggplot(data.frame(density = d_ppc), aes(x = density)) +
  #   geom_histogram(bins = 30, color = "black", fill = "lightblue") +
  #   geom_vline(xintercept = ss$d_obs, color = "red", size = 1) +
  #   labs(
  #     title = paste0(pct_label, "% PPC density"),
  #     x = expression(d == e / (K[D] * K[G]))
  #   ) +
  #   theme_minimal()
  # ggsave(here(plots_dir, paste0("ppc_density_", pct_label, "pct.pdf")), ppc_hist)

  # Bayes factor for sigma_A < delta
  post0 <- mean(draws_df$sigma_A < delta)
  post1 <- mean(draws_df$sigma_A > delta)
  prior0 <- pbeta(delta, 0.05, 1)
  prior1 <- 1 - prior0
  bf <- (post0 / post1) * (prior1 / prior0)
  log10_bf <- log10(bf)

  # store summaries for later plotting
  sigma_summaries <- bind_rows(
    sigma_summaries,
    tibble(
      pct        = 100 * p,
      sigma_A    = draws_df$sigma_A,
      sigma_B    = draws_df$sigma_B
    )
  )
  bf_summary <- bind_rows(
    bf_summary,
    tibble(
      pct     = 100 * p,
      BF      = bf,
      log10BF = log10_bf
    )
  )
}

# ---- Plot posterior densities of sigma_A and sigma_B across sample sizes ----
# we’ll facet by parameter for clarity
sigma_long <- sigma_summaries %>%
  pivot_longer(c(sigma_A, sigma_B), names_to = "parameter", values_to = "value")

dens_all <- ggplot(sigma_long, aes(x = value, color = factor(pct))) +
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  labs(
    color = "% data",
    title = "Posterior densities across subsample sizes"
  ) +
  theme_minimal()
ggsave(here(plots_dir, "sigma_densities_all.pdf"), dens_all, width = 8, height = 4)

# ---- Plot Bayes factor vs. percentage ----
bf_plot <- ggplot(bf_summary, aes(x = pct, y = log10BF)) +
  geom_line() +
  geom_point() +
  labs(
    x = "% data",
    y = "log10 Bayes Factor",
    title = bquote("Bayes factor for " ~ sigma[A] < .(delta))
  ) +
  theme_minimal()
ggsave(here(plots_dir, "bayes_factor_vs_pct.pdf"), bf_plot, width = 6, height = 4)
