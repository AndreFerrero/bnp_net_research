# R script: Stan sampling for Pitman–Yor model with visualization and summary

# Load libraries
library(rstan)
library(bayesplot)
library(ggplot2)
library(here)

# Set options
options(
  mc.cores = 5,
  rstan.auto_write = TRUE
)

source("code/functions/1_py_sample_lab.R")

# Simulate data from PYP
set.seed(42)
alpha_true <- 5
sigma_true <- 0.7

res <- py_sample_lab(
  N            = 1e3,
  alpha        = alpha_true,
  sigma        = sigma_true,
  base_sampler = function() rnorm(1)
)

#---------------------------------------------
# 1) Compile Stan model (once)
#---------------------------------------------
stan_path <- here("code/estimation", "py_sampler.stan")
stan_mod  <- stan_model(file = stan_path)

#---------------------------------------------
# 2) Prepare data and run sampling
#---------------------------------------------

gamma_moments <- function(shape, rate) {
  mean <- shape / rate
  variance <- shape / (rate^2)
  sd <- sqrt(variance)
  
  list(mean = mean, variance = variance, sd = sd)
}

beta_moments <- function(a, b) {
  mean <- a / (a + b)
  variance <- (a * b) / ((a + b)^2 * (a + b + 1))
  sd <- sqrt(variance)
  list(mean = mean, variance = variance, sd = sd)
}

inv_logit <- function(x) {
  p = 1/(1+exp(-x))
  return(p)
}

gamma_moments(3, 0.4)
beta_moments(0.5,5)
inv_logit(0)

stan_data <- list(
  K       = res$n_clusters,
  n_j     = res$counts,
  prior_alpha = c(3, 0.4),
  prior_sigma = c(1, 1)
)

# Fit model
fit <- sampling(
  object = stan_mod,
  data   = stan_data,
  chains = 4,
  iter   = 8000,
  warmup = 1000,
  seed   = 42,
  thin   = 1
)

#---------------------------------------------
# 3) Extract posterior samples
#---------------------------------------------
post <- as.data.frame(rstan::extract(fit, pars = c("alpha", "sigma")))

#---------------------------------------------
# Create plot saving directory
#---------------------------------------------
est_pics_path <- here("res", "pics", "estimation")
if (!dir.exists(est_pics_path)) dir.create(est_pics_path, recursive = TRUE)

#---------------------------------------------
# 4) Diagnostics: traceplots, densities, ACF
#---------------------------------------------
color_scheme_set("brightblue")

# Traceplots
traceplot <- mcmc_trace(as.array(fit), pars = c("alpha", "sigma")) +
  ggtitle("Traceplots") +
  theme(legend.position = "top")
# ggsave(filename = file.path(est_pics_path, "traceplot_alpha_sigma.png"),
#        plot = traceplot, width = 6, height = 5, dpi = 300, bg = "white")

# Density plots
post_alpha <- mcmc_dens_overlay(fit, pars = "alpha") +
  ggtitle("Posterior Density of alpha") +
  geom_vline(xintercept = alpha_true, color = "darkred",
             linetype = "dashed", linewidth = 0.8) +
  theme(legend.position = "top")

post_sigma <- mcmc_dens_overlay(fit, pars = "sigma") +
  ggtitle("Posterior Density of sigma") +
  geom_vline(xintercept = sigma_true, color = "darkorange2",
             linetype = "dashed", linewidth = 0.8) +
  theme(legend.position = "top")

densities_plot <- gridExtra::grid.arrange(post_alpha, post_sigma, ncol = 2)
# ggsave(filename = file.path(est_pics_path, "density_plots_alpha_sigma.png"),
#        plot = densities_plot, width = 6, height = 5, dpi = 300, bg = "white")

# ACF plot
acf_plot <- mcmc_acf_bar(as.array(fit), pars = c("alpha", "sigma")) +
  ggtitle("ACF")
# ggsave(filename = file.path(est_pics_path, "acf_plot_alpha_sigma.png"),
#        plot = acf_plot, width = 6, height = 5, dpi = 300, bg = "white")

# Joint posterior plot
joint_plot <- ggplot(post, aes(x = alpha, y = sigma)) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Joint posterior samples", x = "alpha", y = "sigma")
# ggsave(filename = file.path(est_pics_path, "joint_posterior_samples.png"),
#        plot = joint_plot, width = 6, height = 5, dpi = 300, bg = "white")


post_summary <- summary(fit, pars = c("alpha", "sigma"),
                            probs = c(0.025, 0.5, 0.975))$summary
print(post_summary)