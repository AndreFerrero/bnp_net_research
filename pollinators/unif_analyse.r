library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(tidyverse)
library(posterior)

poll_dir <- here("pollinators/mpl025")

# IMPORT DATA
# Read the CSV file
# (Change the file path to where you've saved your file)
data <- read.csv(here(poll_dir, "M_PL_025.csv"),
  check.names = FALSE, row.names = 1
)

# Convert the adj matrix to an edge list
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

(d_obs <- nrow(edges) / (nrow(plant) * nrow(poll)))


# ESTIMATION

# options(
#   mc.cores = parallel::detectCores(),
#   rstan.auto_write = TRUE
# )

# stan_folder <- here("code", "stan")
# unif_ppc_path <- here(stan_folder, "unif_ppc.stan")
# unif_ppc_mod <- stan_model(file = unif_ppc_path)

# unif_ppc_data <- list(
#   K_A = nrow(plant),
#   K_B = nrow(poll),
#   n_A = plant$counts,
#   n_B = poll$counts,
#   prior_alpha_A = c(3, 0.4),
#   prior_alpha_B = c(3, 0.4),
#   eps = 0.1,
#   e_obs = sum(edges$weight)
# )

# # Fit model
# unif_ppc_fit <- sampling(
#   object = unif_ppc_mod,
#   data   = unif_ppc_data,
#   chains = 4,
#   iter   = 20000,
#   warmup = 1000,
#   seed   = 42,
#   thin   = 10,
#   control = list(adapt_delta = 0.999)
# )

# save(unif_ppc_fit, file = here(poll_dir, "full_unif_ppc_fit.Rdata"))

# LOAD MODEL FIT

load(here(poll_dir, "wo_repl", "unif_fit", "unif_ppc_fit_100pct.Rdata"))

# ——————————————————————————————————————————————————————————————
# 4. Extract and analyze posterior‐predictive density
# ——————————————————————————————————————————————————————————————
draws <- as_draws_df(fit)
d_ppc <- draws$density_ppc

# Summaries of the PPC distribution
d_quantile <- quantile(d_ppc, probs = c(0.025, 0.5, 0.975))
print(d_quantile)
cat("Observed density:", round(d_obs, 4), "\n")

# Posterior p-value
sum(d_ppc >= d_obs) / length(d_ppc)


# Plot the PPC histogram with observed line
ggplot(data.frame(density = d_ppc), aes(x = density)) +
  geom_histogram(bins = 30, color = "black", fill = "lightblue") +
  geom_vline(xintercept = d_obs, color = "red", size = 1) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

ggsave(
  filename = here(poll_dir, "wo_repl", "unif_plots", "full_ppc.pdf"),
  width = 5,
  height = 4
)
