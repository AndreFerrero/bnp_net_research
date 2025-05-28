# this script goes in the cluster
library(tidyverse)
library(here)
library(rstan)

poll_dir <- here("poll")

data <- read.csv(here(poll_dir, "M_PL_006.csv"),
  check.names = FALSE, row.names = 1
)

full_edges <- data |>
  rownames_to_column(var = "plant") |>
  pivot_longer(
    cols = -plant,
    names_to = "pollinator",
    values_to = "weight"
  ) |>
  filter(weight > 0)

total_weight <- sum(full_edges$weight)
cat("Total edge‐weight in network:", total_weight, "\n")

sub_repl <- function(edges, subn) {
  # Sample indices with replacement, weighted by original weights
  sampled_indices <- sample(seq_len(nrow(edges)), size = subn, replace = TRUE, prob = edges$weight)

  # Count the number of times each row (edge) is sampled
  edge_counts <- table(sampled_indices)

  # Reconstruct the sub_edges with updated weights
  sub_edges <- edges[as.integer(names(edge_counts)), ]
  sub_edges$weight <- as.integer(edge_counts)

  # Summarise by plant and pollinator
  sub_plant <- sub_edges |>
    dplyr::group_by(plant) |>
    dplyr::summarise(counts = sum(weight), .groups = "drop")

  sub_poll <- sub_edges |>
    dplyr::group_by(pollinator) |>
    dplyr::summarise(counts = sum(weight), .groups = "drop")

  list(
    sub_plant = sub_plant,
    sub_poll  = sub_poll,
    sub_edges = sub_edges,
    e_obs     = subn
  )
}


# make use of all cores
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Stan model
stan_folder <- here("stan")
beta_ppc_mod <- stan_model(file = here(stan_folder, "beta_ppc.stan"))

# define subn values
reps <- c(10^4, 10^5, 10^6)

# loop over each percentage
for (subn in reps) {
  cat(sprintf("subn = %d\n", subn))

  sub_edges <- sub_repl(full_edges, subn)

  stan_data <- list(
    K_A           = nrow(sub_edges$sub_plant),
    K_B           = nrow(sub_edges$sub_poll),
    n_A           = sub_edges$sub_plant$counts,
    n_B           = sub_edges$sub_poll$counts,
    prior_alpha_A = c(3, 0.4),
    prior_alpha_B = c(3, 0.4),
    eps           = 0.1,
    e_obs         = sub_edges$e_obs
  )

  fit <- sampling(
    object  = beta_ppc_mod,
    data    = stan_data,
    chains  = 4,
    iter    = 20000,
    warmup  = 1000,
    thin    = 10,
    seed    = 42,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )

  check_hmc_diagnostics(fit)

  subn_label <- sprintf("%d", subn)
  out_file <- here(
    "poll", "fits",
    paste0("beta_01_1_ppc_fit_", subn_label, ".Rdata")
  )

  save(fit, file = out_file)
}
