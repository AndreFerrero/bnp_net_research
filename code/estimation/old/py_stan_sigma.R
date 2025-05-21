# R script: Stan sampling for Pitman–Yor model with visualization and summary

# Load libraries
library(rstan)
library(bayesplot)
library(ggplot2)
library(here)

# Set options
options(
  mc.cores = 5,
  rstan.auto_write = TRUE
)

source("code/functions/1_py_sample_lab.R")

# Simulate data from PYP
set.seed(42)
res <- py_sample_lab(
  N = 1e6,
  alpha = 5,
  sigma = 0.3,
  base_sampler = function() rnorm(1)
)

#---------------------------------------------
# 1) Compile Stan model (once)
#---------------------------------------------
# stan_path <- here("code/estimation", "py_sampler_sigma.stan")
# stan_mod  <- stan_model(file = stan_path)

stan_path_eta <- here("code/estimation", "py_sampler_eta.stan")
stan_mod_eta  <- stan_model(file = stan_path_eta)

#---------------------------------------------
# 2) Prepare data and run sampling
#---------------------------------------------

beta_moments <- function(a, b) {
  mean <- a / (a + b)
  variance <- (a * b) / ((a + b)^2 * (a + b + 1))
  list(mean = mean, variance = variance)
}

beta_moments(1,1)

inv_logit <- function(x) {
  p = 1/(1+exp(-x))
  return(p)
}

# stan_data <- list(
#   K       = res$n_clusters,
#   n_j     = res$counts,
#   a_sigma = 0.8,
#   b_sigma = 1.2,
#   alpha_fixed = 5
# )

stan_data_eta <- list(
  K       = res$n_clusters,
  n_j     = res$counts,
  mu = 0,
  tau = 1,
  alpha_fixed = 5
)

# # Fit model
# fit <- sampling(
#   object = stan_mod,
#   data   = stan_data,
#   chains = 4,
#   iter   = 2000,
#   warmup = 1000,
#   seed   = 42
# )

fit_eta <- sampling(
  object = stan_mod_eta,
  data   = stan_data_eta,
  chains = 4,
  iter   = 8000,
  warmup = 1000,
  seed   = 42
)
#---------------------------------------------
# 3) Extract posterior samples
#---------------------------------------------
# post <- as.data.frame(rstan::extract(fit, pars = c("sigma")))

#---------------------------------------------
# 4) Diagnostics: traceplots, densities, ACF
#---------------------------------------------
# color_scheme_set("brightblue")
# 
# # Traceplots
mcmc_trace(fit_eta, pars = c("sigma"))
# 
# Density plots
mcmc_dens_overlay(fit_eta, pars = c("sigma"))
# 
# ACF plots
mcmc_acf(fit_eta, pars = c("sigma"))

#---------------------------------------------
# 5) Posterior summaries and intervals
#---------------------------------------------
# post_summary <- summary(fit, pars = c("sigma"),
#                         probs = c(0.025, 0.5, 0.975))$summary
# print(post_summary)

post_summary_eta <- summary(fit_eta, pars = c("sigma", "eta"),
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_summary_eta)


# MAP estimation

fit_map <- optimizing(
  object    = stan_mod_eta,
  data      = stan_data_eta,
  as_vector = FALSE
)

# Extract the MAP on eta and back‐transform
eta_map   <- fit_map$par$eta
sigma_map <- 1 / (1 + exp(-eta_map))
cat("MAP estimate σ =", sigma_map, "\n")
