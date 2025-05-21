source('code/funs/py_sample.R')
source("code/funs/sample_net.R")

# Load libraries
library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(grid)

# Set options
options(
  mc.cores = 5,
  rstan.auto_write = TRUE
)

set.seed(42)

alpha_true = c(5,5)
sigma_true = c(0.7, 0.7)

net <- sample_net(1e4,
                  alpha = alpha_true,
                  sigma = sigma_true)

#---------------------------------------------
# 1) Compile Stan model (once)
#---------------------------------------------
stan_folder <- here("code", "stan")
unif_path <- here(stan_folder, "unif_model.stan")
unif_mod  <- stan_model(file = unif_path)

stan_data <- list(
  K_A = net$xA$K,
  K_B = net$xB$K,
  n_A = net$xA$active_counts,
  n_B = net$xB$active_counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  prior_sigma_A = c(1, 1),
  prior_sigma_B = c(1, 1)
)

# Fit model
unif_fit <- sampling(
  object = unif_mod,
  data   = stan_data,
  chains = 4,
  iter   = 2000,
  warmup = 1000,
  seed   = 42,
  thin   = 1
)


est_folder <- here("code","estimation")

save(unif_fit, file = here(est_folder, "unif_fit.Rdata"))

#---------------------------------------------
# 3) Extract posterior samples
#---------------------------------------------
post <- as.data.frame(rstan::extract(fit, pars = c("alpha_A",
                                                   "alpha_B",
                                                   "sigma_A",
                                                   "sigma_B")))

#---------------------------------------------
# Create plot saving directory
#---------------------------------------------
est_pics_path <- here("res", "pics", "estimation", "py_net")
if (!dir.exists(est_pics_path)) dir.create(est_pics_path, recursive = TRUE)

#---------------------------------------------
# 4) Diagnostics: traceplots, densities, ACF
#---------------------------------------------
color_scheme_set("brightblue")

# Traceplots
traceplot <- mcmc_trace(as.array(fit), pars = c("alpha_A",
                                                "alpha_B",
                                                "sigma_A",
                                                "sigma_B")) +
  ggtitle("Traceplots") +
  theme(legend.position = "top")
# ggsave(filename = file.path(est_pics_path, "trace_50py_net0707.png"),
#        plot = traceplot, width = 7, height = 6, dpi = 300, bg = "white")

# Density plots
post_alphaA <- mcmc_dens_overlay(as.array(fit), pars = c("alpha_A")) +
  geom_vline(xintercept = alpha_true[1], color = "darkred",
             linetype = "dashed", linewidth = 0.8) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(x = expression(alpha[A]))

post_alphaB <- mcmc_dens_overlay(as.array(fit), pars = c("alpha_B")) +
  geom_vline(xintercept = alpha_true[2], color = "darkred",
             linetype = "dashed", linewidth = 0.8) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(x = expression(alpha[B]))

post_sigmaA <- mcmc_dens_overlay(as.array(fit), pars = "sigma_A") +
  geom_vline(xintercept = sigma_true[1], 
             color = "darkorange2", linetype = "dashed", linewidth = 0.8) +
  labs(x = expression(sigma[A])) +
  theme(legend.position = "none")

post_sigmaB <- mcmc_dens_overlay(as.array(fit), pars = "sigma_B") +
  geom_vline(xintercept = sigma_true[2], 
             color = "darkorange2", linetype = "dashed", linewidth = 0.8) +
  labs(x = expression(sigma[B])) +
  theme(legend.position = "none")

post_alpha <- gridExtra::grid.arrange(post_alphaA, post_alphaB,
                                      ncol = 2)
post_sigma <- gridExtra::grid.arrange(post_sigmaA, post_sigmaB,
                                      ncol = 2)

main_title <- textGrob("Posterior densities", gp = gpar(fontsize = 14), 
                       hjust = 0.5, vjust = 1)  # Centered title

# Arrange the plots and add the main title
densities_plot <- gridExtra::grid.arrange(main_title, post_alpha, post_sigma,
                                          nrow = 3, heights = c(0.1, 1, 1))
# ggsave(filename = file.path(est_pics_path, "post_50py_net0707.png"),
#        plot = densities_plot, width = 7, height = 6, dpi = 300, bg = "white")

# ACF plot
acf_plot <- mcmc_acf_bar(as.array(fit), pars = c("alpha_A",
                                                 "alpha_B",
                                                 "sigma_A",
                                                 "sigma_B")) +
  ggtitle("ACF")
# ggsave(filename = file.path(est_pics_path, "acf_50py_net0707.png"),
#        plot = acf_plot, width = 7, height = 6, dpi = 300, bg = "white")

# Joint posterior plot
# joint_plot <- ggplot(post, aes(x = alpha, y = sigma)) +
#   geom_point(alpha = 0.3) +
#   theme_minimal() +
#   labs(title = "Joint posterior samples", x = "alpha", y = "sigma")
# ggsave(filename = file.path(est_pics_path, "joint_posterior_samples.png"),
#        plot = joint_plot, width = 7, height = 6, dpi = 300, bg = "white")


post_summary <- summary(fit, pars = c("alpha_A",
                                      "alpha_B",
                                      "sigma_A",
                                      "sigma_B"),
                        probs = c(0.025, 0.5, 0.975))$summary
print(post_summary)