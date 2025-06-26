library(rstan)
library(here)
library(posterior)
library(bayesplot)
library(ggplot2)
library(tidyverse)

# pics folder
unif_0_02_pics_folder <- here("res", "pics", "estimation", "0_02_unif_ppc")
unif_0_07_pics_folder <- here("res", "pics", "estimation", "0_07_unif_ppc")

# est folder
est_folder <- here("code", "estimation")

load(here(est_folder, "0_02_unif_ppc_fit.Rdata"))
unif_0_02_ppc_fit <- unif_ppc_fit
load(here(est_folder, "0_07_unif_ppc_fit.Rdata"))
unif_0_07_ppc_fit <- unif_ppc_fit
rm(unif_ppc_fit)

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
# unif_ppc_path <- here(stan_folder, "unif_ppc.stan")
# unif_ppc_mod <- stan_model(file = unif_ppc_path)

# # Data for Stan
# unif_ppc_data <- list(
#   K_A = net$xA$K,
#   K_B = net$xB$K,
#   n_A = net$xA$active_counts,
#   n_B = net$xB$active_counts,
#   prior_alpha_A = c(3, 0.4),
#   prior_alpha_B = c(3, 0.4),
#   prior_sigma_A = c(1, 1),
#   prior_sigma_B = c(1, 1),
#   e_obs = nrow(net$edges)
# )

# # Fit the model
# unif_ppc_fit <- sampling(
#   object = unif_ppc_mod,
#   data = unif_ppc_data,
#   chains = 4,
#   iter = 4000,
#   warmup = 1000,
#   seed = 42,
#   thin = 2,
#   control = list(adapt_delta = 0.95)
# )

# Save the fit
# save(unif_ppc_fit, file = here(est_folder, "0_07_unif_ppc_fit.Rdata"))

# check_hmc_diagnostics(unif_ppc_fit)

# Posterior summary
unif_summary <- summary(unif_ppc_fit, probs = c(0.025, 0.5, 0.975))$summary
print(round(unif_summary, digits = 3))

bayes_plots <- function(fit_object, color_set, net, pics_folder, save = FALSE) {
  color_scheme_set(color_set)

  if (color_set == "blue") {
    hist_color <- "#b3cde0"
  } else {
    hist_color <- "#DCBCBC"
  }

  param_labels <- c(
    alpha_A = "alpha[A]",
    alpha_B = "alpha[B]",
    sigma_A = "sigma[A]",
    sigma_B = "sigma[B]"
  )

  # ACF plot
  (acf_plot <- mcmc_acf_bar(fit_object,
    pars = names(param_labels),
    facet_args = list(
      labeller = labeller(
        .variables = NULL,
        .default = label_parsed,
        .cols = param_labels
      )
    )
  ) +
    theme(legend.position = "none") +
    hline_at(0.2, linetype = 2, size = 0.15, color = "gray") +
    scale_y_continuous(breaks = c(0.2)))

  # Trace plot
  (trace_plot <- mcmc_trace(fit_object,
    pars = names(param_labels),
    facet_args = list(
      nrow = 2,
      labeller = labeller(
        .variables = NULL, .default = label_parsed,
        .cols = param_labels
      )
    )
  ) +
    theme(legend.position = "none")
  )

  # Posterior densities
  (dens_plot <- mcmc_dens_overlay(fit_object,
    pars = names(param_labels),
    facet_args = list(
      nrow = 2,
      labeller = labeller(
        .variables = NULL, .default = label_parsed,
        .cols = param_labels
      )
    )
  ) +
    theme(legend.position = "none")
  )

  # PPC observed network density

  # Convert to posterior draws
  unif_ppc_draws <- as_draws_df(fit_object)

  d_obs <- nrow(unique(net$edges)) / (net$xA$K * net$xB$K)
  unif_d_ppc <- unif_ppc_draws$density_ppc

  hist_data <- data.frame(density = unif_d_ppc) %>%
    mutate(bin = cut(density, breaks = 30)) %>%
    group_by(bin) %>%
    summarise(
      count = n(),
      x = mean(range(density[bin == bin[1]])) # bin center
    ) %>%
    mutate(
      color = ifelse(x >= d_obs, hist_color, "grey89") # conditional color
    )

  # Build the plot
  ppc_plot <- ggplot(hist_data, aes(x = x, y = count, fill = color)) +
    geom_col(color = "grey22", width = diff(range(unif_d_ppc)) / 30) +
    scale_fill_identity() + # use actual hex values from 'fill'
    annotate(
      "segment",
      x = d_obs, xend = d_obs,
      y = 0, yend = Inf,
      color = color_set,
      size = 0.5,
      linetype = "dashed"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_blank())

  if (save) {
    ggsave(
      filename = here(pics_folder, "unif_ppc_acf.pdf"),
      plot = acf_plot,
      device = "pdf",
      width = 6,
      height = 4
    )

    ggsave(
      filename = here(pics_folder, "unif_ppc_trace.pdf"),
      plot = trace_plot,
      device = "pdf",
      width = 6,
      height = 4
    )

    ggsave(
      filename = here(pics_folder, "unif_ppc_post.pdf"),
      plot = dens_plot,
      width = 4,
      height = 3
    )

    ggsave(
      filename = here(pics_folder, "unif_ppc_dens.pdf"),
      plot = ppc_plot,
      width = 4,
      height = 3
    )
  }
}

bayes_plots(unif_0_02_ppc_fit, "red", net_0_02, unif_0_02_pics_folder, save = T)
bayes_plots(unif_0_07_ppc_fit, "blue", net_0_07, unif_0_07_pics_folder, save = T)



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
