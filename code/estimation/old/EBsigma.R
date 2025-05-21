
## Sampling ####
source("code/functions/1_py_sample_lab.R")

alpha_fixed = 5
sigma_true = 0.3

# Simulate data from PYP
set.seed(42)
res <- py_sample_lab(
  N = 1e4,
  alpha = alpha_fixed,
  sigma = sigma_true,
  base_sampler = function() rnorm(1)
)

# Rising Pochhammer: (a)_m = Gamma(a+m) / Gamma(a)
log_pochhammer <- function(a, m) {
  lgamma(a + m) - lgamma(a)
}

## EPPF with sigma ####

# Pitman–Yor EPPF log-likelihood
eppf_lp <- function(sigma, n_j, alpha) {
  K <- length(n_j)
  N <- sum(n_j)
  
  # 1) sum_{i=1}^{K-1} log(alpha + i * sigma)
  term1 <- sum(log(alpha + (1:(K-1)) * sigma))
  
  # 2) - log_pochhammer(alpha + 1, N - 1)
  term2 <- - log_pochhammer(alpha + 1, N - 1)
  
  # 4) sum_j log_pochhammer(1 - sigma, n_j[j] - 1)
  term3 <- sum(sapply(n_j, function(nj) log_pochhammer(1 - sigma, nj - 1)))
  
  term1 + term2 + term3
}

#––– Create a grid of sigma values in (0,1) –––
n_grid <- 1000
sigma_grid <- seq(1e-4, 1 - 1e-4, length.out = n_grid)

loglik_vals <- sapply(sigma_grid,
                      function(s) eppf_lp(s, res$counts, alpha_fixed))

#––– 4) Plot with ggplot2 –––
library(ggplot2)
df_grid <- data.frame(sigma = sigma_grid, loglik = loglik_vals)

ggplot(df_grid, aes(x = sigma, y = loglik)) +
  geom_line(size = 1) +
  labs(
    x     = expression(sigma),
    y     = "log EPPF",
    title = "Pitman–Yor EPPF log-likelihood over σ"
  ) +
  theme_minimal()

## Logit transformation ####

eppf_lp_transf <- function(x, n_j, alpha) {
  K <- length(n_j)
  N <- sum(n_j)
  
  # inverse transformation
  sigma <- 1 / (1 + exp(-x))
  
  # 1) sum_{i=1}^{K-1} log(alpha + i * sigma)
  term1 <- sum(log(alpha + (1:(K-1)) * sigma))
  
  # 2) - log_pochhammer(alpha + 1, N - 1)
  term2 <- - log_pochhammer(alpha + 1, N - 1)
  
  # 3) sum_j log_pochhammer(1 - sigma, n_j[j] - 1)
  term3 <- sum(sapply(n_j, function(nj) log_pochhammer(1 - sigma, nj - 1)))
  
  term1 + term2 + term3
}

# grid over x
x_grid <- seq(-5, 5, length.out = 1000)

# compute log-EPPF at each x
loglik_x <- sapply(x_grid, eppf_lp_transf, n_j = res$counts,
                   alpha = alpha_fixed)

df <- data.frame(x = x_grid, loglik = loglik_x)

ggplot(df, aes(x = x, y = loglik)) +
  geom_line(size = 1) +
  labs(
    x     = expression(x),
    y     = "log EPPF",
    title = "Pitman–Yor log-EPPF over x"
  ) +
  theme_minimal()

## Empirical Bayes with optim() ####

# Optimization
neg_eppf <- function(sigma, n_j, alpha) {
  # negates the log-likelihood
  -eppf_lp(sigma, n_j, alpha)
}

neg_eppf_transf <- function(x, n_j, alpha_fixed) {
  -eppf_lp_transf(x, n_j, alpha_fixed)
}


opt_res <- optim(
  par     = 0.5,                # initial guess
  fn      = neg_eppf,          
  gr      = NULL,               # let optim approximate gradient
  lower   = 1e-8,               # enforce sigma > 0
  upper   = 1 - 1e-8,           # enforce sigma < 1
  method  = "L-BFGS-B",
  control = list(fnscale = 1),  # fnscale=1 since we're already providing neg-log-lik
  n_j     = res$counts,                # extra data passed to neg_eppf
  alpha   = alpha_fixed
)


sigma_hat  <- opt_res$par

cat("Empirical‐Bayes σ̂ =", sigma_hat, "\n",
    "True sigma =", sigma_true, "\n")


opt_res_transf <- optim(
  par    = 0,               # init x = 0 → σ = 0.5
  fn     = neg_eppf_transf,
  method = "BFGS",
  n_j    = res$counts,
  alpha_fixed = alpha_fixed,
  control = list(fnscale = 1)
)

sigma_hatx <- 1 / (1 + exp(-opt_res_transf$par))

cat("Empirical‐Bayes σ̂ with transformation =", sigma_hatx, "\n",
    "True sigma =", sigma_true, "\n")
