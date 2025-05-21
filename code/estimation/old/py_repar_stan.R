library(rstan)
library(ggplot2)
library(coda)
library(dplyr)
library(tidyr)
library(here)

source("functions/1_py_sample_lab.R")
rstan_options(auto_write = TRUE)
options(mc.cores = 5)

# Simulate data
set.seed(42)
res <- py_sample_lab(
  N = 1e5,
  alpha = 5,
  sigma = 0.6,
  base_sampler = function() rnorm(1)
)

n_j <- as.integer(res$counts)

# Prepare data list
stan_data <- list(
  K = res$n_clusters,
  n_j = n_j
)

# Compile Stan model
stan_path <- here("code/estimation", "py_sampler_reparam.stan")
model <- stan_model(stan_path)

# Fit model
fit <- sampling(
  model,
  data = stan_data,
  chains = 4,
  iter = 20000,
  warmup = 1000,
  seed = 42
)

#---------------------------------------------
# 3) Extract posterior samples
#---------------------------------------------
post <- as.data.frame(rstan::extract(fit, pars = c("alpha", "sigma")))

#---------------------------------------------
# 4) Diagnostics: traceplots, densities, ACF
#---------------------------------------------
color_scheme_set("brightblue")

# Traceplots
mcmc_trace(as.array(fit), pars = c("alpha", "sigma")) +
  ggtitle("Traceplots for alpha and sigma")

# Density plots
mcmc_dens_overlay(as.array(fit), pars = c("alpha", "sigma")) +
  ggtitle("Posterior Densities: alpha and sigma")

# ACF plots
mcmc_acf(as.array(fit), pars = c("alpha", "sigma")) +
  ggtitle("ACF for alpha and sigma")

#---------------------------------------------
# 5) Posterior summaries and intervals
#---------------------------------------------
post_summary <- summary(fit, pars = c("alpha", "sigma"), probs = c(0.025, 0.5, 0.975))$summary
print(post_summary)

# Convert to tibble for nicer formatting
# Use tibble::rownames_to_column or load tibble
library(tibble)
post_sum_df <- as.data.frame(post_summary) %>%
  rownames_to_column(var = "parameter") %>%
  dplyr::select(parameter, mean, `2.5%`, `50%`, `97.5%`, n_eff, Rhat)

print(post_sum_df)

#---------------------------------------------
ggplot(post, aes(x = alpha, y = sigma)) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Joint posterior of alpha and sigma",
       x = "alpha", y = "sigma")
