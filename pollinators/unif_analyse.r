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


(total_weight <- sum(edges$weight))
(binary_weight <- nrow(edges))


# compare with theoretical network
source("code/funs/py_sample.R")
source("code/funs/sample_net.R")

# True parameters
alpha_true <- c(5, 5)
sigma_true <- c(0, 0.1)

tnet <- sample_net(total_weight, alpha_true, sigma_true)

nrow(unique(tnet$edges))
nrow(tnet$edges) == total_weight

nrow(unique(tnet$edges)) / (tnet$xA$K * tnet$xB$K)

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
  width = 7,
  height = 4
)


# Parameters to extract
params <- c("alpha_A", "alpha_B", "sigma_A", "sigma_B")

# Extract diagnostic summary
extract_diagnostics <- function(fit, params) {
  summ <- summary(fit, pars = params)$summary
  data.frame(
    Rhat = round(summ[, "Rhat"], 3),
    n_eff = round(summ[, "n_eff"], 0)
  )
}

diag <- extract_diagnostics(fit, params)

# Get number of divergences
get_divergences <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
}

div <- get_divergences(fit)

# Combine into one data frame
diag_df <- data.frame(
  Parameter = params,
  Rhat = diag$Rhat,
  Neff = diag$n_eff
)

# Add math formatting
diag_df$Parameter <- recode(diag_df$Parameter,
  "alpha_A" = "\\alpha_A",
  "alpha_B" = "\\alpha_B",
  "sigma_A" = "\\sigma_A",
  "sigma_B" = "\\sigma_B"
)

# Print LaTeX table
cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Convergence diagnostics for the estimation on real data, including $\\widehat{R}$, effective sample size $n_{\\mathrm{eff}}$, and number of divergent transitions.}\n")
cat("\\label{tab:real_diag}\n")
cat("\\begin{tabular}{lcc}\n")
cat("\\toprule\n")
cat("Parameter & $\\widehat{R}$ (1) & $n_{\\mathrm{eff}}$ (1) & $\\widehat{R}$ (2) & $n_{\\mathrm{eff}}$ (2) \\\\\n")
cat("\\midrule\n")

for (i in 1:nrow(diag_df)) {
  cat(paste0(
    "$", diag_df$Parameter[i], "$ & ",
    diag_df$Rhat[i], " & ", diag_df$Neff[i]" \\\\\n"
  ))
}

# Add divergences row
cat("\\midrule\n")
cat(paste0("Divergences & \\multicolumn{2}{c}{"div"}\\\\\n"))

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
