# This script simulates one single full network and then
# estimates the basic model on subnets of different sizes.

library(rstan)
library(here)

# Load helper functions (expects sample_net to exist and return at least $edges)
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

# Desired subnet edge sizes (prefixes)
edge_sizes <- c(1e2, 1e3, 1e4, 1e5) # full is max

# Stan models
stan_folder <- here("code", "stan")
mod_path <- here(stan_folder, "basic_model.stan")
ppc_mod_path <- here(stan_folder, "basic_ppc.stan")

mod <- stan_model(file = mod_path)
ppc_mod <- stan_model(file = ppc_mod_path)

# Estimation folder
est_folder <- here("code", "estimation")
if (!dir.exists(est_folder)) dir.create(est_folder, recursive = TRUE)

# Prepare data from an edges data.frame
prepare_data_from_edges <- function(edges_df, include_e_obs = FALSE) {
  # Compute frequencies of each unique node on sides A and B
  freq_A <- table(edges_df$X_A)
  freq_B <- table(edges_df$X_B)

  data_list <- list(
    K_A = length(freq_A),
    K_B = length(freq_B),
    n_A = as.integer(freq_A),
    n_B = as.integer(freq_B),
    prior_alpha_A = c(3, 0.4),
    prior_alpha_B = c(3, 0.4),
    prior_sigma_A = c(1, 1),
    prior_sigma_B = c(1, 1)
  )
  if (include_e_obs) {
    data_list$e_obs <- nrow(edges_df)
  }
  return(data_list)
}

# Toggle PPC
ppc <- FALSE # set to TRUE to run PPC on the 1e4-edge subnet

# === Main simulation & estimation workflow ===
for (sigma_label in names(sigma_values)) {
  sigma <- sigma_values[[sigma_label]]

  # Simulate only the largest network once
  n_full <- max(edge_sizes)
  cat(sprintf("Simulating FULL network for n = %g, sigma = %s\n", n_full, sigma_label))
  full_net <- sample_net(n_full, alpha = alpha_true, sigma = sigma)
  edges_full <- full_net$edges # expect data.frame with X_A and X_B

  # Loop over prefix subnet sizes
  for (n_edges in edge_sizes) {
    # Take prefix of edges
    edges_sub <- head(edges_full, n_edges)

    # Prepare Stan data (basic model)
    stan_data_basic <- prepare_data_from_edges(edges_sub, include_e_obs = FALSE)

    cat(sprintf("Fitting BASIC model for subnet n = %g, sigma = %s\n", n_edges, sigma_label))
    fit <- sampling(
      object = mod,
      data = stan_data_basic,
      chains = 4,
      iter = 5000,
      warmup = 1000,
      seed = 42,
      thin = 1
    )

    save_file <- here(est_folder, sprintf("fit_subnet_n_%g_sigma_%s.Rdata", n_edges, sigma_label))
    save(fit, file = save_file)
    rm(fit)

    # PPC for the 1e4 subnet if enabled
    if (n_edges == 1e4 && ppc == TRUE) {
      cat(sprintf("Running PPC model for subnet n = %g, sigma = %s\n", n_edges, sigma_label))
      stan_data_ppc <- prepare_data_from_edges(edges_sub, include_e_obs = TRUE)

      fit_ppc <- sampling(
        object = ppc_mod,
        data = stan_data_ppc,
        chains = 4,
        iter = 5000,
        warmup = 1000,
        seed = 42,
        thin = 2
      )

      ppc_save_file <- here(est_folder, sprintf("fit_ppc_subnet_n_%g_sigma_%s.Rdata", n_edges, sigma_label))
      save(fit_ppc, file = ppc_save_file)
      rm(fit_ppc)
    }
  }
}
