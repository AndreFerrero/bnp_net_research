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

sub_wo_repl <- function(edges, subn) {
  set.seed(1234)

  # original counts
  plant <- edges |>
    group_by(plant) |>
    summarise(counts = sum(weight))
  poll <- edges |>
    group_by(pollinator) |>
    summarise(counts = sum(weight))

  # explode into tickets and sample exactly subn of them
  tickets <- rep(seq_len(nrow(edges)), times = edges$weight)
  sel <- sample(tickets, size = subn, replace = FALSE)

  sub_weights_df <- as.data.frame(table(sel), stringsAsFactors = FALSE) |>
    transmute(
      row_index = as.integer(sel),
      weight    = as.integer(Freq)
    )

  sub_edges <- edges |>
    slice(sub_weights_df$row_index) |>
    mutate(weight = sub_weights_df$weight)

  sub_plant <- sub_edges |>
    group_by(plant) |>
    summarise(counts = sum(weight))

  sub_poll <- sub_edges |>
    group_by(pollinator) |>
    summarise(counts = sum(weight))

  list(
    sub_plant = sub_plant,
    sub_poll  = sub_poll,
    sub_edges = sub_edges,
    e_obs     = sum(sub_edges$weight)
  )
}

# make use of all cores
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Stan model
stan_folder <- here("stan")
beta_ppc_mod <- stan_model(file = here(stan_folder, "beta_ppc.stan"))

# define the fractions of the data you want to try
percentages <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

# loop over each percentage
for (p in percentages) {
  subn <- round(total_weight * p)
  cat(sprintf("→ %3.0f%% of total weight → subn = %d\n", p * 100, subn))

  sub_edges <- sub(full_edges, subn)

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

  pct_label <- sprintf("%03d", round(p * 100))
  out_file <- here(
    "poll", "fits",
    paste0("beta_05_1_ppc_fit_", pct_label, "pct.Rdata")
  )
  save(fit, file = out_file)
}
