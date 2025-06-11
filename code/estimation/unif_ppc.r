library(rstan)
library(here)
library(posterior)
library(bayesplot)
library(ggplot2)

# pics folder
unif_ppc_pics_folder <- here("res", "pics", "estimation", "0_07_unif_ppc")

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
sigma_true <- c(0, 0.7)

# Simulate network data
net <- sample_net(1e4,
  alpha = alpha_true,
  sigma = sigma_true
)

# Compile Stan model
unif_ppc_path <- here(stan_folder, "unif_ppc.stan")
unif_ppc_mod <- stan_model(file = unif_ppc_path)

# Data for Stan
unif_ppc_data <- list(
  K_A = net$xA$K,
  K_B = net$xB$K,
  n_A = net$xA$active_counts,
  n_B = net$xB$active_counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  prior_sigma_A = c(1, 1),
  prior_sigma_B = c(1, 1),
  e_obs = nrow(net$edges)
)

# Fit the model
unif_ppc_fit <- sampling(
  object = unif_ppc_mod,
  data = unif_ppc_data,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 42,
  thin = 2,
  control = list(adapt_delta = 0.95)
)

# Save the fit
save(unif_ppc_fit, file = here(est_folder, "0_07_unif_ppc_fit.Rdata"))

check_hmc_diagnostics(unif_ppc_fit)

# Posterior summary
unif_summary <- summary(unif_ppc_fit, probs = c(0.025, 0.5, 0.975))$summary
print(round(unif_summary, digits = 3))

# Trace plot
trace_plot <- mcmc_trace(unif_ppc_fit, pars = c("alpha_A", "alpha_B", "sigma_A", "sigma_B")) +
  theme(legend.position = "top")

ggsave(
  filename = here(unif_ppc_pics_folder, "unif_ppc_trace.pdf"),
  plot = trace_plot,
  device = "pdf",
  width = 6,
  height = 4
)

# ACF plot
acf_plot <- mcmc_acf_bar(unif_ppc_fit, pars = c("alpha_A", "alpha_B", "sigma_A", "sigma_B"))

ggsave(
  filename = here(unif_ppc_pics_folder, "unif_ppc_acf.pdf"),
  plot = acf_plot,
  device = "pdf",
  width = 6,
  height = 4
)

# Posterior densities
dens_plot <- mcmc_dens_overlay(unif_ppc_fit, pars = c("alpha_A", "alpha_B", "sigma_A", "sigma_B")) +
  theme(legend.position = "none")

ggsave(
  filename = here(unif_ppc_pics_folder, "unif_ppc_post.pdf"),
  plot = dens_plot,
  width = 4,
  height = 3
)

# PPC observed network density

# Convert to posterior draws
unif_ppc_draws <- as_draws_df(unif_ppc_fit)

d_obs <- nrow(unique(net$edges)) / (net$xA$K * net$xB$K)
unif_d_ppc <- unif_ppc_draws$density_ppc

# Summary of PPC
d_quantile <- quantile(unif_d_ppc, probs = c(0.025, 0.5, 0.975))
print(d_quantile)
cat("Observed density:", round(d_obs, 4), "\n")

# PPC histogram
ppc_plot <- ggplot(data.frame(density = unif_d_ppc), aes(x = density)) +
  geom_histogram(bins = 30, color = "black", fill = "lightblue") +
  geom_vline(xintercept = d_obs, color = "red", size = 1) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank())

ggsave(
  filename = here(unif_ppc_pics_folder, "unif_ppc_dens.pdf"),
  plot = ppc_plot,
  width = 4,
  height = 3
)

# Bayes Factor
delta <- 0.05
post_0 <- mean(unif_ppc_draws$sigma_A < delta)
post_1 <- mean(unif_ppc_draws$sigma_A > delta)
prior_0 <- pbeta(delta, 1, 1)
prior_1 <- pbeta(delta, 1, 1, lower.tail = FALSE)

bf <- post_0 / post_1 * prior_1 / prior_0
log10_bf <- log10(bf)

cat("Bayes Factor (BF10):", round(bf, 3), "\n")
cat("log10(BF10):", round(log10_bf, 3), "\n")
