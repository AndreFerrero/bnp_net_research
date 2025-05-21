source('code/funs/py_sample.R')
source("code/funs/sample_net.R")

# Load libraries
library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(grid)

beta_moments <- function(a, b) {
  mean <- a / (a + b)
  variance <- (a * b) / ((a + b)^2 * (a + b + 1))
  list(mean = mean, variance = variance)
}

beta_moments(0.05,1)

# Set options
options(
  mc.cores = 4,
  rstan.auto_write = TRUE
)

set.seed(42)

alpha_true = c(5,5)
sigma_true = c(0, 0.7)

net <- sample_net(1e4,
                  alpha = alpha_true,
                  sigma = sigma_true)

stan_folder <- here("code", "stan")
spike_path <- here(stan_folder, "spike_model.stan")

spike_mod <- stan_model(file = spike_path)

spike_data <- list(
  K_A = net$xA$K,
  K_B = net$xB$K,
  n_A = net$xA$active_counts,
  n_B = net$xB$active_counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  eps = 0.05
)

# Fit model
spike_fit <- sampling(
  object = spike_mod,
  data   = spike_data,
  chains = 4,
  iter   = 4000,
  warmup = 1000,
  seed   = 42,
  thin   = 1,
  control = list(adapt_delta = 0.99)
)

est_folder <- here("code","estimation")

save(spike_ppc_fit, file = here(est_folder, "spike_fit_eps_005.Rdata"))

post_summary <- summary(spike_fit,
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_summary |> round(digits = 3))

spike_folder <- here("res", "pics", "estimation", "spike")

mcmc_trace(spike_fit, pars = c("alpha_A",
                               "alpha_B",
                               "sigma_A",
                               "sigma_B")) +
  ggtitle("Spike-slab prior on sigma") +
  theme(legend.position = "top")

mcmc_trace(spike_fit, pars = c("alpha_A",
                               "alpha_B",
                               "sigma_A",
                               "sigma_B",
                               "w_A",
                               "w_B")) +
  ggtitle("Spike-slab prior on sigma") +
  theme(legend.position = "top")

ggsave(here(spike_folder, "trace.pdf"),
       width = 12,
       height = 10,
       bg = "white",
       dpi = 600)

mcmc_acf_bar(spike_fit, pars = c("alpha_A",
                               "alpha_B",
                               "sigma_A",
                               "sigma_B",
                               "w_A",
                               "w_B")) +
  ggtitle("ACF") +
  theme(legend.position = "top")

ggsave(here(spike_folder, "acf.pdf"),
       width = 12,
       height = 10,
       bg = "white",
       dpi = 600)

mcmc_pairs(spike_fit, pars = c("sigma_A", "w_A"))
mcmc_pairs(spike_fit, pars = c("sigma_B", "w_B"))

mcmc_dens_overlay(spike_fit, pars = c("alpha_A",
                                      "alpha_B",
                                      "sigma_A",
                                      "sigma_B",
                                      "w_A",
                                      "w_B"))

ggsave(here(spike_folder, "post.pdf"),
       width = 12,
       height = 10,
       bg = "white",
       dpi = 600)

## Region of Practical Equivalence test (ROPE) ####
draws <- as.data.frame(spike_fit)

delta <- 0.05

# P(sigma_A < delta | Y)
sum(draws$sigma_A < delta)/nrow(draws)

## E[w|Y] ####

# extract posterior means for w_A and w_B
w_means <- post_summary[c("w_A", "w_B"), "mean"]

# compute pseudo-Bayes factors
pseudo_BF <- w_means / (1 - w_means)

# print results
print(round(pseudo_BF, 3))
