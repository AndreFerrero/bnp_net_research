library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(tidyverse)
library(posterior)
library(bipartite)

poll_dir <- here("pollinators/mpl025")
wo_repl_dir <- here(poll_dir, "wo_repl")
plots_dir <- here(wo_repl_dir, "unif_plots")
fits_dir <- here(wo_repl_dir, "unif_fit")

# Import adjacency matrix
data <- read.csv(here(poll_dir, "M_PL_025.csv"),
  check.names = FALSE, row.names = 1
)

# Edge list and observed d
edges <- data |>
  rownames_to_column(var = "plant") |>
  pivot_longer(
    cols = -plant,
    names_to = "pollinator",
    values_to = "weight"
  ) |>
  filter(weight > 0)

plant <- edges |>
  group_by(plant) |>
  summarise(counts = sum(weight))
poll <- edges |>
  group_by(pollinator) |>
  summarise(counts = sum(weight))

d_obs <- nrow(edges) / (nrow(plant) * nrow(poll))
d_obs

web <- ifelse(data > 0, 1, 0)
density <- function(web) {
  # web: binary adjacency matrix (plants rows x pollinators cols)
  total_links <- sum(web)
  possible_links <- nrow(web) * ncol(web)
  d <- total_links / possible_links
  return(d)
}

set.seed(123)
n_sims <- 1000

# EE null model: equiprobable
d_EE <- replicate(n_sims, {
  web_sim <- matrix(rbinom(length(web), size = 1, prob = d_obs),
    nrow = nrow(web)
  )
  sum(web_sim) / (nrow(plant) * nrow(poll))
})

load(here(fits_dir, "unif_ppc_fit_100pct.Rdata"))
draws_df <- as_draws_df(fit)

# posterior predictive draws for d
d_ppc <- draws_df$density_ppc

d_df <- tibble(
  d = c(d_EE, d_ppc),
  model = c(
    rep("EE null", length(d_EE)),
    rep("Bayesian PPC", length(d_ppc))
  )
)

ggplot(d_df, aes(x = d, fill = model)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = d_obs, colour = "black", linetype = "dashed", size = 1) +
  labs(
    title = "Observed d vs Null Model vs Bayesian PPC",
    subtitle = "Dashed line = observed d",
    x = "d (connectance)",
    y = "Density"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "EE null" = "#F8766D",
    "Bayesian PPC" = "#7CAE00"
  ))

ppp_bayes <- mean(d_ppc >= d_obs)
ppp_EE <- mean(d_EE >= d_obs)

ppp_bayes
ppp_EE
