source("functions/1_py_sample_lab.R")

set.seed(42)
res <- py_sample_lab(
  N            = 1e3,
  alpha        = 20,
  sigma        = 0.3,
  base_sampler = function() rnorm(1)
)

# Extract cluster sizes
cluster_sizes <- res$counts

# Define log-eppf with fixed alpha ----
log_eppf <- function(par, N, alpha_fixed) {
  x <- par[1]
  sigma <- 1 / (1 + exp(-x))
  alpha <- alpha_fixed
  
  if (alpha <= -sigma) return(-Inf)
  
  K <- length(N)
  n <- sum(N)
  
  term1 <- sum(log(alpha + sigma * (0:(K - 1))))
  term2 <- sum(log(alpha + 1 + (0:(n - 2))))
  term3 <- sum(lgamma(N - sigma) - lgamma(1 - sigma))
  
  term1 - term2 + term3
}

obj_fn <- function(par, N, alpha_fixed) -log_eppf(par, N, alpha_fixed)

# Optimize ----
alpha_fixed <- 20  # You can change this
init_par <- 0  # logit of initial sigma (0.5)
res_opt <- optim(
  par     = init_par,
  fn      = obj_fn,
  N       = cluster_sizes,
  alpha_fixed = alpha_fixed,
  method  = "BFGS",
  control = list(fnscale = 1)
)

# Back-transform ----
x_hat     <- res_opt$par
sigma_hat <- 1 / (1 + exp(-x_hat))

cat(sprintf("Empirical Bayes estimate:\n  sigma = %.4f\n", sigma_hat))
