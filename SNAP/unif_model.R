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
unif_path <- here(stan_folder, "unif_model.stan")
unif_mod <- stan_model(file = unif_path)

unif_data <- list(
  K_A = K_D,
  K_B = K_G,
  n_A = n_D,
  n_B = n_G,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  prior_sigma_A = c(1,1),
  prior_sigma_B = c(1,1)
)

unif_fit <-  sampling(
  object = unif_mod,
  data   = unif_data,
  chains = 4,
  iter   = 6000,
  warmup = 2000,
  seed   = 42,
)

check_hmc_diagnostics(unif_fit)

save(unif_fit, file = here(snap_folder, "unif_fit.Rdata"))

post_unif_summary <- summary(unif_fit,
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_unif_summary |> round(digits = 3))

mcmc_trace(unif_fit, pars = c("alpha_A",
                               "alpha_B",
                               "sigma_A",
                               "sigma_B")) +
  ggtitle("Uniform prior on sigma") +
  theme(legend.position = "top")

ggsave(filename = here(snap_folder, "pics","unif_trace.png"),
       width = 7, height = 6, bg = "white")

mcmc_acf_bar(unif_fit, pars = c("alpha_A",
                                 "alpha_B",
                                 "sigma_A",
                                 "sigma_B")) +
  ggtitle("Uniform prior on sigma") +
  theme(legend.position = "top")

ggsave(filename = here(snap_folder, "pics","unif_acf.png"),
       width = 7, height = 6, bg = "white")

mcmc_dens_overlay(unif_fit,
                  pars = c("alpha_A",
                           "alpha_B",
                           "sigma_A",
                           "sigma_B")) +
  theme(legend.position = "top") +
  ggtitle("Uniform prior on sigma")

ggsave(filename = here(snap_folder, "pics","unif_post.png"),
       width = 7, height = 6, bg = "white")