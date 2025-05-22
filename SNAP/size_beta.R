#!/usr/bin/env Rscript

library(rstan)
library(bayesplot)
library(ggplot2)
library(posterior)
library(readr)
library(here)
library(dplyr)
library(tidyr)

options(mc.cores = parallel::detectCores())

#——— Paths ——————————————————————————————————————————————————————
snap_folder    <- here("SNAP")
stan_folder    <- here("stan")
results_folder <- here(snap_folder, "results")
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

#——— Data & Model —————————————————————————————————————————————
data_path    <- here(snap_folder, "gene_desease.gz")
edges        <- read_tsv(data_path, col_names = TRUE)
beta_ppc_mod <- stan_model(file = here(stan_folder, "beta_ppc.stan"))
eps <- 0.1

#——— Sample sizes ——————————————————————————————————————————————
subn_values <- c(100, 500, 1000, 5000, 10000)
# deltas for Bayes factor
delta_values <- c(0.001, 0.005, 0.01)

#——— Storage for combined plots & BF table ————————————————————————————————————
all_draws <- vector("list", length(subn_values))
names(all_draws) <- paste0("n", subn_values)
# Create a BF table with all combinations of subn and delta
bf_table <- expand.grid(subn = subn_values, delta = delta_values) %>%
  arrange(subn, delta) %>%
  mutate(bayes_factor = NA_real_)

#——— Loop ——————————————————————————————————————————————————————
for (i in seq_along(subn_values)) {
  subn <- subn_values[i]
  message("▶ subn = ", subn)
  run_dir <- file.path(results_folder, paste0("subn_", subn))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1) Sample data
  set.seed(123 + subn)
  edges_iid <- edges[sample(nrow(edges), subn), ]
  colnames(edges_iid) <- c("Disease", "Gene")
  
  K_D  <- length(unique(edges_iid$Disease))
  K_G  <- length(unique(edges_iid$Gene))
  n_D  <- table(edges_iid$Disease)
  n_G  <- table(edges_iid$Gene)
  e    <- nrow(edges_iid)
  d_obs <- e / (K_D * K_G)
  
  stan_data <- list(
    K_A           = K_D,
    K_B           = K_G,
    n_A           = n_D,
    n_B           = n_G,
    prior_alpha_A = c(3, 0.4),
    prior_alpha_B = c(3, 0.4),
    eps           = eps,
    e_obs         = e
  )
  
  # 2) Fit model (parallel chains)
  fit <- sampling(
    object  = beta_ppc_mod,
    data    = stan_data,
    chains  = 4, iter = 6000, warmup = 2000, thin = 2,
    seed    = 42 + subn,
    control = list(adapt_delta = 0.999)
  )
  save(fit, file = file.path(run_dir, "beta_ppc_fit.Rdata"))
  
  # 3) Extract draws
  draws <- as_draws_df(fit)
  d_ppc <- draws$density_ppc
  q     <- quantile(d_ppc, probs = c(0.025, 0.5, 0.975))
  
  # 4) Bayes factors for each delta
  for (j in seq_along(delta_values)) {
    delta <- delta_values[j]
    post_0  <- mean(draws$sigma_A < delta)
    post_1  <- mean(draws$sigma_A > delta)
    prior_0 <- pbeta(delta, eps, 1)
    prior_1 <- pbeta(delta, eps, 1, lower.tail = FALSE)
    bf      <- (post_0 / post_1) * (prior_1 / prior_0)
    idx     <- which(bf_table$subn == subn & bf_table$delta == delta)
    bf_table$bayes_factor[idx] <- bf
  }
  
  # 5) Save summary
  summary_lines <- c(
    paste0("Observed density      : ", round(d_obs, 4)),
    paste0("2.5% / 50% / 97.5%   : ", round(q[1],4), " / ", round(q[2],4), " / ", round(q[3],4)),
    ""
  )
  for (delta in delta_values) {
    bf <- bf_table %>% filter(subn == subn, delta == delta) %>% pull(bayes_factor)
    summary_lines <- c(summary_lines,
                       paste0("BF (σ_A < ", delta, ") = ", signif(bf, 3))
    )
  }
  writeLines(summary_lines, con = file.path(run_dir, "density_and_bf_summary.txt"))
  
  # 6) Individual plots
  p_trace <- mcmc_trace(fit, pars = c("sigma_A", "sigma_B", "density_ppc")) +
    ggtitle(paste("Trace – subn =", subn)) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))
  ggsave(file.path(run_dir, "trace_plot.pdf"), p_trace, width = 10, height = 6, bg = "white")
  
  p_ppc <- ggplot(data.frame(d = d_ppc), aes(x = d)) +
    geom_histogram(bins = 30, color = "black", fill = "lightblue") +
    geom_vline(xintercept = d_obs, color = "red", size = 1) +
    labs(
      title = paste("PPC – subn =", subn),
      x = expression(paste("Density = e / (K[D]*K[G])")),
      y = "Count"
    ) + theme_minimal() +
    theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))
  ggsave(file.path(run_dir, "ppc_density_plot.pdf"), p_ppc, width = 8, height = 5, bg = "white")
  
  # 7) Store draws
  all_draws[[i]] <- draws %>%
    select(sigma_A, sigma_B, density_ppc) %>%
    mutate(subn = subn) %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")
  
  rm(fit, draws, d_ppc)
  gc()
}

#——— Combined plots —————————————————————————————————————————————————
combo_df <- bind_rows(all_draws)

# A) σ_A
p_sigmaA <- combo_df %>%
  filter(parameter == "sigma_A") %>%
  ggplot(aes(x = value, color = factor(subn))) +
  geom_density() +
  labs(color = "subn", title = "Posterior density of σ_A", x = "σ_A", y = "Density") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))
ggsave(file.path(results_folder, "combined_sigmaA.pdf"), p_sigmaA, width = 8, height = 5, bg = "white")

# B) σ_B
p_sigmaB <- combo_df %>%
  filter(parameter == "sigma_B") %>%
  ggplot(aes(x = value, color = factor(subn))) +
  geom_density() +
  labs(color = "subn", title = "Posterior density of σ_B", x = "σ_B", y = "Density") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))
ggsave(file.path(results_folder, "combined_sigmaB.pdf"), p_sigmaB, width = 8, height = 5, bg = "white")

# C) density_ppc
p_ppc_all <- combo_df %>%
  filter(parameter == "density_ppc") %>%
  ggplot(aes(x = value, color = factor(subn))) +
  geom_density() +
  labs(color = "subn", title = "Posterior-predictive density",
       x = expression(paste("Density = e/(K[D]*K[G])")),
       y = "Density") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
ggsave(file.path(results_folder, "combined_ppc_density.pdf"), p_ppc_all,
       width = 8, height = 5, bg = "white")

# D) Bayes factor vs. subn
p_bf <- ggplot(bf_table, aes(x = subn, y = bayes_factor, group = factor(delta), color = factor(delta))) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = subn_values) +
  labs(
    title = "Bayes Factor (σ_A < δ) vs. sample size",
    x = "subn",
    y = "Bayes factor",
    color = expression(delta)
  ) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))
ggsave(file.path(results_folder, "bf_vs_subn.pdf"), p_bf, width = 8, height = 5, bg = "white")

message("✅ All runs complete. Results under: ", results_folder)
