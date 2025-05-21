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
beta_path <- here(stan_folder, "beta_model.stan")

beta_mod <- stan_model(file = beta_path)

beta_data <- list(
  K_A = net$xA$K,
  K_B = net$xB$K,
  n_A = net$xA$active_counts,
  n_B = net$xB$active_counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  eps = 0.05
)

# Fit model
beta_fit <- sampling(
  object = beta_mod,
  data   = beta_data,
  chains = 4,
  iter   = 4000,
  warmup = 1000,
  seed   = 42,
  thin   = 2,
  control = list(adapt_delta = 0.999)
)

est_folder <- here("code","estimation")

save(beta_fit, file = here(est_folder, "beta_fit_eps_005.Rdata"))

post_summary <- summary(beta_fit,
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_summary |> round(digits = 3))

beta_folder <- here("res", "pics", "estimation", "beta")

mcmc_trace(beta_fit, pars = c("alpha_A",
                               "alpha_B",
                               "sigma_A",
                               "sigma_B")) +
  ggtitle("beta-slab prior on sigma") +
  theme(legend.position = "top")

ggsave(here(spike_folder, "trace.pdf"),
       width = 12,
       height = 10,
       bg = "white",
       dpi = 600)

mcmc_acf_bar(beta_fit, pars = c("alpha_A",
                                 "alpha_B",
                                 "sigma_A",
                                 "sigma_B")) +
  ggtitle("ACF") +
  theme(legend.position = "top")

ggsave(here(spike_folder, "acf.pdf"),
       width = 12,
       height = 10,
       bg = "white",
       dpi = 600)

mcmc_pairs(beta_fit, pars = c("sigma_A", "w_A"))
mcmc_pairs(beta_fit, pars = c("sigma_B", "w_B"))

mcmc_dens_overlay(beta_fit, pars = c("alpha_A",
                                      "alpha_B",
                                      "sigma_A",
                                      "sigma_B"))

ggsave(here(spike_folder, "post.pdf"),
       width = 12,
       height = 10,
       bg = "white",
       dpi = 600)

## BAYES FACTOR ####
draws <- as.data.frame(beta_fit)

delta <- 0.05

# P(sigma_A < delta | Y)
post_0 <- mean(draws$sigma_A < delta)
post_1 <- mean(draws$sigma_A > delta)

prior_0 <- pbeta(delta, 0.05, 1)
prior_1 <- pbeta(delta, 0.05, 1, lower.tail = F)

(bf <- post_0/post_1 * prior_1/prior_0)
