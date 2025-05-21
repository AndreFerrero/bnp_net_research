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
spike_ppc_path <- here(stan_folder, "spike_ppc.stan")
spike_ppc_mod <- stan_model(file = spike_ppc_path)

spike_ppc_data <- list(
  K_A = K_D,
  K_B = K_G,
  n_A = n_D,
  n_B = n_G,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  eps = 0.1,
  e_obs = e
)

# Fit model
spike_ppc_fit <- sampling(
  object = spike_ppc_mod,
  data   = spike_ppc_data,
  chains = 4,
  iter   = 6000,
  warmup = 2000,
  seed   = 42,
  thin   = 2,
  control = list(adapt_delta = 0.999)
)


save(spike_ppc_fit, here(snap_folder, "spike_ppc_fit.Rdata"))

mcmc_trace(spike_ppc_fit, pars = c("alpha_A",
                               "alpha_B",
                               "sigma_A",
                               "sigma_B",
                               "w_A",
                               "w_B",
                               "density_ppc")) +
  ggtitle("Spike-slab prior on sigma") +
  theme(legend.position = "top")

# ——————————————————————————————————————————————————————————————
# 4. Extract and analyze posterior‐predictive density
# ——————————————————————————————————————————————————————————————

library(posterior)

draws <- as_draws_df(spike_ppc_fit)
d_ppc <- draws$density_ppc

# 4.3 Summaries of the PPC distribution
d_quantile <- quantile(d_ppc, probs = c(0.025, 0.5, 0.975))
print(d_quantile)
cat("Observed density:", round(d_obs, 4), "\n")

# 4.4 Plot the PPC histogram with observed line
ggplot(data.frame(density = d_ppc), aes(x = density)) +
  geom_histogram(bins = 30, color = "black", fill = "lightblue") +
  geom_vline(xintercept = d_obs, color = "red", size = 1) +
  labs(
    title = "PPC - Networks Density",
    x = expression(paste("Density = ", e / (K[D] * K[G]))),
    y = "Count"
  ) +
  theme_minimal()

rm(spike_ppc_fit)
