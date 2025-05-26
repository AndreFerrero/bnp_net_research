library(rstan)
library(bayesplot)
library(ggplot2)
library(here)

poll_dir <- here("poll")
source(here(poll_dir, "import.r"))

options(
  mc.cores = parallel::detectCores(),
  rstan.auto_write = TRUE
)

stan_folder <- here("stan")
beta_ppc_path <- here(stan_folder, "beta_ppc.stan")
beta_ppc_mod <- stan_model(file = beta_ppc_path)

beta_ppc_data <- list(
  K_A = nrow(plant),
  K_B = nrow(poll),
  n_A = plant$counts,
  n_B = poll$counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  eps = 0.1,
  e_obs = sum(edges$weight)
)

# Fit model
beta_ppc_fit <- sampling(
  object = beta_ppc_mod,
  data   = beta_ppc_data,
  chains = 4,
  iter   = 20000,
  warmup = 1000,
  seed   = 42,
  thin   = 10,
  control = list(adapt_delta = 0.999)
)

save(beta_ppc_fit, file = here(poll_dir, "full_beta_ppc_fit.Rdata"))

# load(here(poll_dir, "full_beta_ppc_fit.Rdata"))

post_summary <- summary(beta_ppc_fit,
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_summary |> round(digits = 3))

mcmc_trace(beta_ppc_fit, pars = c("alpha_A",
                                  "alpha_B",
                                  "sigma_A",
                                  "sigma_B")) +
  ggtitle("Beta prior on sigma") +
  theme(legend.position = "top")

ggsave(filename = here("pollinators", "mpl006", "full_trace_thin_10.pdf"))

mcmc_acf_bar(beta_ppc_fit, pars = c("alpha_A",
                                  "alpha_B",
                                  "sigma_A",
                                  "sigma_B")) +
  ggtitle("Beta prior on sigma") +
  theme(legend.position = "top")

ggsave(filename = here("pollinators", "mpl006", "full_acf_thin_10.pdf"))

mcmc_dens_overlay(beta_ppc_fit, pars = c("alpha_A",
                                  "alpha_B",
                                  "sigma_A",
                                  "sigma_B")) +
  ggtitle("Beta prior on sigma") +
  theme(legend.position = "top")

ggsave(filename = here("pollinators", "mpl006", "full_dens_thin_10.pdf"))

# ——————————————————————————————————————————————————————————————
# 4. Extract and analyze posterior‐predictive density
# ——————————————————————————————————————————————————————————————

library(posterior)
draws <- as_draws_df(beta_ppc_fit)
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
    x = expression(d == e / (K[D] * K[G])),
    y = NULL
  ) +
  theme_minimal()

ggsave(filename = here("pollinators", "mpl006", "full_ppc_thin_10.pdf"))

## BAYES FACTOR ####
delta <- 0.05

# P(sigma_A < delta | Y)
(post_0 <- mean(draws$sigma_A < delta))
post_1 <- mean(draws$sigma_A > delta)

prior_0 <- pbeta(delta, 0.05, 1)
prior_1 <- pbeta(delta, 0.05, 1, lower.tail = FALSE)

(bf <- post_0 / post_1 * prior_1 / prior_0)
log10(bf)
