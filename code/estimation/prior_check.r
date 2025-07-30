set.seed(123)

source("code/funs/py_sample.R")
source("code/funs/sample_net.R")

# Define prior sampling functions
sample_prior_alpha <- function() {
  rgamma(2, shape = 3, rate = 0.4)
}

# Two versions for sigma
sample_prior_sigma_beta <- function() {
  rbeta(2, shape1 = 1.5, shape2 = 3)
}

sample_prior_sigma_unif <- function() {
  runif(2, min = 0, max = 1) # uniform prior
}

n_prior_sims <- 500

# --- Beta(2,3) prior ---
prior_samples_beta <- replicate(n_prior_sims, {
  alpha <- sample_prior_alpha()
  sigma <- sample_prior_sigma_beta()

  net <- sample_net(N = 2225, alpha = alpha, sigma = sigma)

  nrow(unique(net$edges)) / (net$xA$K * net$xB$K)
})

# --- Uniform prior ---
prior_samples_unif <- replicate(n_prior_sims, {
  alpha <- sample_prior_alpha()
  sigma <- sample_prior_sigma_unif()

  net <- sample_net(N = 2000, alpha = alpha, sigma = sigma)

  nrow(unique(net$edges)) / (net$xA$K * net$xB$K)
})

library(tidyverse)

prior_df <- tibble(
  density = c(prior_samples_unif, prior_samples_beta),
  sigma_prior = c(
    rep("Uniform(0,1)", n_prior_sims),
    rep("Beta(1.5,3)", n_prior_sims)
  )
)

ggplot(prior_df, aes(x = density, fill = sigma_prior)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Prior Predictive Distribution of Density",
    x = "Density (d)", y = "Density (pdf)"
  ) +
  theme_minimal()
