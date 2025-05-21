
library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
snap_folder = here("SNAP")

source(here(snap_folder, "import.R"))

options(
  mc.cores = 4,
  rstan.auto_write = TRUE
)

stan_folder <- here("code","stan")
spike_path <- here(stan_folder, "spike_model.stan")
spike_mod <- stan_model(file = spike_path)

spike_data <- list(
  K_A = K_D,
  K_B = K_G,
  n_A = n_D,
  n_B = n_G,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  eps = 0.1
)

# Fit model
spike_fit <- sampling(
  object = spike_mod,
  data   = spike_data,
  chains = 4,
  iter   = 6000,
  warmup = 2000,
  seed   = 42,
  thin   = 2,
  control = list(adapt_delta = 0.999)
)


check_hmc_diagnostics(spike_fit)

post_summary <- summary(spike_fit,
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_summary |> round(digits = 3))

mcmc_trace(spike_fit, pars = c("alpha_A",
                               "alpha_B",
                               "sigma_A",
                               "sigma_B",
                               "w_A",
                               "w_B")) +
  ggtitle("Spike-slab prior on sigma") +
  theme(legend.position = "top")

ggsave(filename = here(snap_folder, "pics","spike_trace.png"),
       width = 7, height = 6, bg = "white")

mcmc_acf_bar(spike_fit, pars = c("alpha_A",
                                 "alpha_B",
                                 "sigma_A",
                                 "sigma_B",
                                 "w_A",
                                 "w_B")) +
  ggtitle("Spike-slab prior on sigma") +
  theme(legend.position = "top")

ggsave(filename = here(snap_folder, "pics","spike_acf.png"),
       width = 7, height = 6, bg = "white")

mcmc_dens_overlay(spike_fit,
                  pars = c("alpha_A",
                           "alpha_B",
                           "sigma_A",
                           "sigma_B",
                           "w_A",
                           "w_B")) +
  theme(legend.position = "top") +
  ggtitle("Spike-slab prior on sigma")

ggsave(filename = here(snap_folder, "pics","spike_post.png"),
       width = 7, height = 6, bg = "white")

mcmc_dens_overlay(spike_fit, pars = "sigma_B") +
  xlim(0, 2e-05)

mcmc_pairs(spike_fit, pars = c("sigma_A", "w_A"))
mcmc_pairs(spike_fit, pars = c("sigma_B", "w_B"))

## Region of Practical Equivalence test (ROPE) ####
draws <- as.data.frame(spike_fit)

delta <- 0.01

# P(sigma_A < delta | Y)
sum(draws$sigma_A < delta)/nrow(draws)

## E[w|Y] ####

# extract posterior means for w_A and w_B
w_means <- post_summary[c("w_A", "w_B"), "mean"]

# compute pseudo-Bayes factors
pseudo_BF <- w_means / (1 - w_means)

# print results
print(round(pseudo_BF, 3))
