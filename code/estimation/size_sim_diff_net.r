# This script simulates networks of different sizes
# and estimates the basic model on them.

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

set.seed(1234)

# True parameters
alpha_true <- c(5, 5)
sigma_values <- list(
  "0_02" = c(0, 0.2),
  "0_07" = c(0, 0.7)
)

# Edge sizes to simulate
edge_sizes <- c(1e2, 1e3, 1e4, 1e5)

# Stan models
stan_folder <- here("code", "stan")
mod_path <- here(stan_folder, "basic_model.stan")
ppc_mod_path <- here(stan_folder, "basic_ppc.stan")

mod <- stan_model(file = mod_path)
ppc_mod <- stan_model(file = ppc_mod_path)

# Estimation folder
est_folder <- here("code", "estimation")

# Function to prepare data for Stan
prepare_data <- function(net, include_e_obs = FALSE) {
  data_list <- list(
    K_A = net$xA$K,
    K_B = net$xB$K,
    n_A = net$xA$active_counts,
    n_B = net$xB$active_counts,
    prior_alpha_A = c(3, 0.4),
    prior_alpha_B = c(3, 0.4),
    prior_sigma_A = c(1, 1),
    prior_sigma_B = c(1, 1)
  )
  if (include_e_obs) {
    data_list$e_obs <- nrow(net$edges)
  }
  return(data_list)
}

ppc <- FALSE

# Simulation loop
for (n_edges in edge_sizes) {
  for (sigma_label in names(sigma_values)) {
    sigma <- sigma_values[[sigma_label]]

    # Simulate network
    net <- sample_net(n_edges, alpha = alpha_true, sigma = sigma)

    # Prepare data for basic model (no e_obs)
    stan_data_basic <- prepare_data(net, include_e_obs = FALSE)

    # Fit basic model
    cat(sprintf("Fitting BASIC model for n = %g, sigma = %s\n", n_edges, sigma_label))
    fit <- sampling(
      object = mod,
      data = stan_data_basic,
      chains = 4,
      iter = 5000,
      warmup = 1000,
      seed = 42,
      thin = 1
    )

    # Save result
    save_file <- here(est_folder, sprintf("fit_n_%g_sigma_%s.Rdata", n_edges, sigma_label))
    save(fit, file = save_file)
    rm(fit)

    # Run PPC model only for n = 1e4
    if (n_edges == 1e4 & ppc == TRUE) {
      cat(sprintf("Running PPC model for n = %g, sigma = %s\n", n_edges, sigma_label))

      # PPC needs e_obs
      stan_data_ppc <- prepare_data(net, include_e_obs = TRUE)

      fit_ppc <- sampling(
        object = ppc_mod,
        data = stan_data_ppc,
        chains = 4,
        iter = 5000,
        warmup = 1000,
        seed = 42,
        thin = 2
      )

      ppc_save_file <- here(est_folder, sprintf("fit_ppc_n_%g_sigma_%s.Rdata", n_edges, sigma_label))
      save(fit_ppc, file = ppc_save_file)
    }
  }
}
