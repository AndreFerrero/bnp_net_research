library(rstan)
library(here)

# Load helper functions
source("code/funs/py_sample.R")
source("code/funs/sample_net.R")

# Set options
options(
  mc.cores = parallel::detectCores(),
  rstan.auto_write = TRUE
)

set.seed(42)

# True parameters
alpha_true <- c(5, 5)
sigma_true_0_02 <- c(0, 0.2)
sigma_true_0_07 <- c(0, 0.7)

# Simulate network data
net_0_02 <- sample_net(1e4,
  alpha = alpha_true,
  sigma = sigma_true_0_02
)
net_0_07 <- sample_net(1e4,
  alpha = alpha_true,
  sigma = sigma_true_0_07
)

# Compile Stan model
stan_folder <- here("stan")
mod_path <- here(stan_folder, "basic_model.stan")
mod <- stan_model(file = mod_path)
ppc_mod_path <- here(stan_folder, "basic_ppc.stan")
ppc_mod <- stan_model(file = ppc_mod_path)

# Data for Stan
data_0_02_ppc <- list(
  K_A = net_0_02$xA$K,
  K_B = net_0_02$xB$K,
  n_A = net_0_02$xA$active_counts,
  n_B = net_0_02$xB$active_counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  prior_sigma_A = c(1.5, 3),
  prior_sigma_B = c(1.5, 3),
  e_obs = nrow(net_0_02$edges)
)

data_0_07_ppc <- list(
  K_A = net_0_07$xA$K,
  K_B = net_0_07$xB$K,
  n_A = net_0_07$xA$active_counts,
  n_B = net_0_07$xB$active_counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  prior_sigma_A = c(1.5, 3),
  prior_sigma_B = c(1.5, 3),
  e_obs = nrow(net_0_07$edges)
)

# Fit the model
fit_0_02_ppc <- sampling(
  object = ppc_mod,
  data = data_0_02,
  chains = 4,
  iter = 5000,
  warmup = 1000,
  seed = 42,
  thin = 1
)
check_hmc_diagnostics(fit_0_02)


fit_0_07_ppc <- sampling(
  object = ppc_mod,
  data = data_0_07,
  chains = 4,
  iter = 5000,
  warmup = 1000,
  seed = 42,
  thin = 1
)
check_hmc_diagnostics(fit_0_07)

# Save the fit
est_folder <- here("estimation")
save(fit_0_02, file = here(est_folder, "0_02_fit.Rdata"))
save(fit_0_07, file = here(est_folder, "0_07_fit.Rdata"))
