#!/usr/bin/env Rscript

#— Load libraries —————————————————————————————————————————————
library(rstan)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(stringr)

#— Settings —————————————————————————————————————————————————————
results_folder <- here("SNAP", "results")
subn_dirs      <- list.dirs(results_folder, recursive = FALSE, full.names = TRUE)
delta_values   <- c(0.001, 0.005, 0.01)

#— 1) Rebuild all_draws —————————————————————————————————————————————
all_draws <- lapply(subn_dirs, function(dir) {
  # extract subn from folder name
  subn <- as.integer(str_extract(basename(dir), "\\d+"))
  
  # load fit
  fit_file <- file.path(dir, "beta_ppc_fit.Rdata")
  if (!file.exists(fit_file)) stop("Missing fit in ", dir)
  load(fit_file)  # creates object `fit`
  
  # extract draws
  draws <- as_draws_df(fit) %>%
    select(sigma_A, sigma_B, density_ppc) %>%
    mutate(subn = subn) %>%
    pivot_longer(cols = c(sigma_A, sigma_B, density_ppc),
                 names_to = "parameter", values_to = "value")
  
  return(draws)
})

combo_df <- bind_rows(all_draws)

#— 2) Rebuild BF table from summaries —————————————————————————————————
# Each folder has density_and_bf_summary.txt with lines:
#   BF (σ_A < delta) = ...
bf_list <- lapply(subn_dirs, function(dir) {
  subn <- as.integer(str_extract(basename(dir), "\\d+"))
  sum_file <- file.path(dir, "density_and_bf_summary.txt")
  lines    <- readLines(sum_file)
  # extract BF lines
  bf_lines <- grep("^BF \\(σ_A", lines, value = TRUE)
  
  df <- data.frame(
    subn = subn,
    delta = sapply(bf_lines, function(l) {
      as.numeric(str_extract(l, "(?<=< )[^\\)]+"))
    }),
    bayes_factor = sapply(bf_lines, function(l) {
      as.numeric(str_extract(l, "(?<== )[-0-9e\\.]+"))
    }),
    stringsAsFactors = FALSE
  )
  return(df)
})

bf_table <- bind_rows(bf_list)

#— 3) Reproduce and save combined plots ——————————————————————————————

# ensure white background
theme_set(
  theme_minimal() +
    theme(panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))
)

# A) σ_A densities
p_sigmaA <- combo_df %>%
  filter(parameter == "sigma_A") %>%
  ggplot(aes(x = value, color = factor(subn))) +
  geom_density(size = 0.8) +
  labs(title = expression("Posterior density of" * " " * sigma[A]),
       x     = expression(sigma[A]),
       y     = "Density",
       color = "subn") 

ggsave(file.path(results_folder, "post_sigmaA.pdf"),
       p_sigmaA, width = 8, height = 5, bg = "white")

# B) σ_B densities
p_sigmaB <- combo_df %>%
  filter(parameter == "sigma_A") %>%
  ggplot(aes(x = value, color = factor(subn))) +
  geom_density(size = 0.8) +
  labs(title = expression("Posterior density of" * " " * sigma[B]),
       x     = expression(sigma[B]),
       y     = "Density",
       color = "subn")

ggsave(file.path(results_folder, "post_sigmaB.pdf"),
       p_sigmaB, width = 8, height = 5, bg = "white")

# C) density_ppc densities
p_ppc <- combo_df %>%
  filter(parameter == "density_ppc") %>%
  ggplot(aes(x = value, color = factor(subn))) +
  geom_density(size = 0.8) +
  labs(title = "Posterior-predictive density",
       x     = expression(paste("Density = e/(K[D]*K[G])")),
       y     = "Density",
       color = "subn")

ggsave(file.path(results_folder, "reproduced_ppc_density.pdf"),
       p_ppc, width = 8, height = 5, bg = "white")

# D) Bayes Factor vs. subn
p_bf <- ggplot(bf_table, aes(x = subn, y = bayes_factor,
                             group = factor(delta), color = factor(delta))) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_x_log10(breaks = unique(bf_table$subn)) +
  labs(title = expression("Bayes Factor ("*sigma[A] < delta*") vs. sample size"),
       x     = "subn",
       y     = "Bayes factor",
       color = expression(delta))

ggsave(file.path(results_folder, "reproduced_bf_vs_subn.pdf"),
       p_bf, width = 8, height = 5, bg = "white")

message("✅ Reproduction complete. PDFs written to ", results_folder)
