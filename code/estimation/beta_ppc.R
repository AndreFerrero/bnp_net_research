library(rstan)
library(bayesplot)
library(ggplot2)
library(here)

source('code/funs/py_sample.R')
source("code/funs/sample_net.R")

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

stan_folder <- here("code","stan")
beta_ppc_path <- here(stan_folder, "beta_ppc.stan")
beta_ppc_mod <- stan_model(file = beta_ppc_path)

beta_ppc_data <- list(
  K_A = net$xA$K,
  K_B = net$xB$K,
  n_A = net$xA$active_counts,
  n_B = net$xB$active_counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  eps = 0.05,
  e_obs = nrow(net$edges)
)

# Fit model
beta_ppc_fit <- sampling(
  object = beta_ppc_mod,
  data   = beta_ppc_data,
  chains = 4,
  iter   = 4000,
  warmup = 1000,
  seed   = 42,
  thin   = 2,
  control = list(adapt_delta = 0.999)
)

est_folder <- here("code","estimation")

save(beta_ppc_fit, file = here(est_folder, "beta_ppc_fit_eps_005.Rdata"))
post_summary <- summary(beta_ppc_fit,
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_summary |> round(digits = 3))



mcmc_trace(beta_ppc_fit, pars = c("alpha_A",
                                   "alpha_B",
                                   "sigma_A",
                                   "sigma_B",
                                   "density_ppc")) +
  ggtitle("Beta prior on sigma") +
  theme(legend.position = "top")

mcmc_trace(beta_ppc_fit, pars = c("alpha_A",
                                   "alpha_B",
                                   "sigma_A",
                                   "sigma_B")) +
  ggtitle("Beta prior on sigma") +
  theme(legend.position = "top")

# ——————————————————————————————————————————————————————————————
# 4. Extract and analyze posterior‐predictive density
# ——————————————————————————————————————————————————————————————
d_obs <- nrow(unique(net$edges))/(net$xA$K * net$xB$K)

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

## BAYES FACTOR ####
delta <- 0.005

# P(sigma_A < delta | Y)
(post_0 <- mean(draws$sigma_A < delta))
post_1 <- mean(draws$sigma_A > delta)

prior_0 <- pbeta(delta, 0.05, 1)
prior_1 <- pbeta(delta, 0.05, 1, lower.tail = F)

(bf <- post_0/post_1 * prior_1/prior_0)
