# combined_script.R

# ——————— Load dependencies ———————
library(ggplot2)
library(here)
source("code/functions/1_py_sample_lab.R")

# ——————— True params & sample data ———————
set.seed(42)
alpha_true <- 5
sigma_true <- 0.7

res <- py_sample_lab(
  N            = 1e4,
  alpha        = alpha_true,
  sigma        = sigma_true,
  base_sampler = function() rnorm(1)
)
n_j <- res$counts

# ——————— Helper transforms & pochhammer ———————
inv_logit     <- function(x) 1/(1+exp(-x))
logit         <- function(p) log(p/(1-p))
log_pochhammer <- function(a,m) lgamma(a+m) - lgamma(a)

# ——————— Log-EPPF under (x=logit(σ), y=log(α)) ———————
eppf_lp_transf <- function(par, n_j) {
  x     <- par[1]
  y     <- par[2]
  sigma <- inv_logit(x)
  alpha <- exp(y)
  K <- length(n_j);  N <- sum(n_j)
  t1 <- sum(log(alpha + (1:(K-1))*sigma))
  t2 <- -log_pochhammer(alpha+1, N-1)
  t3 <- sum(vapply(n_j, function(nj)
    log_pochhammer(1 - sigma, nj - 1), numeric(1)))
  t1 + t2 + t3
}
obj_fn <- function(par, n_j) -eppf_lp_transf(par, n_j)

# ——————— Empirical‐Bayes optimization ———————
init_par <- c(logit(0.5), log(5))
opt     <- optim(
  par     = init_par,
  fn      = obj_fn,
  method  = "BFGS",
  control = list(fnscale=1),
  n_j     = n_j
)

x_hat     <- opt$par[1]
y_hat     <- opt$par[2]
sigma_hat <- inv_logit(x_hat)
alpha_hat <- exp(y_hat)

cat(sprintf("Estimated:  σ̂ = %.4f   α̂ = %.4f\n", sigma_hat, alpha_hat))
cat(sprintf("True:       σ  = %.4f   α  = %.4f\n", sigma_true, alpha_true))

# ——————— Base (untransformed) log-EPPF ———————
eppf_lp <- function(sigma, n_j, alpha) {
  # rising Pochhammer
  lpoch <- function(a,m) lgamma(a+m) - lgamma(a)
  K <- length(n_j);  N <- sum(n_j)
  t1 <- sum(log(alpha + (1:(K-1))*sigma))
  t2 <- -lpoch(alpha+1, N-1)
  t3 <- sum(vapply(n_j, function(nj)
    lpoch(1 - sigma, nj - 1), numeric(1)))
  t1 + t2 + t3
}

# ——————— Build all four grids ———————
n_grid <- 1000

# σ ∈ (0,1)
sigma_grid <- seq(1e-4, 1-1e-4, length.out=n_grid)
ll_sigma   <- sapply(sigma_grid, function(s)
  eppf_lp(s, n_j, alpha_hat)
)
df_sigma <- data.frame(coord="sigma", x=sigma_grid, loglik=ll_sigma)

# logit(σ) on ℝ
x_grid <- seq(logit(1e-4), logit(1-1e-4), length.out=n_grid)
sigma_x <- inv_logit(x_grid)
ll_logit_sigma <- sapply(sigma_x, function(s)
  eppf_lp(s, n_j, alpha_hat)
)
df_logit_sigma <- data.frame(coord="logit_sigma", x=x_grid, loglik=ll_logit_sigma)

# α > 0
alpha_grid <- seq(1e-4, alpha_hat*1.5, length.out=n_grid)
ll_alpha   <- sapply(alpha_grid, function(a)
  eppf_lp(sigma_hat, n_j, a)
)
df_alpha <- data.frame(coord="alpha", x=alpha_grid, loglik=ll_alpha)

# log(α) on ℝ
y_grid  <- seq(log(1e-4), log(alpha_hat)+2, length.out=n_grid)
alpha_y <- exp(y_grid)
ll_logalpha <- sapply(alpha_y, function(a)
  eppf_lp(sigma_hat, n_j, a)
)
df_logalpha <- data.frame(coord="logalpha", x=y_grid, loglik=ll_logalpha)

# Combine
df_all <- rbind(df_sigma, df_logit_sigma, df_alpha, df_logalpha)

# ——————— Plot with EB estimates ———————
# Add factor with desired facet order
df_all$coord <- factor(df_all$coord,
                       levels = c("sigma", "logit_sigma", "alpha", "logalpha"))

# EB estimates for vertical lines
vlines_df <- data.frame(
  coord = factor(c("sigma", "logit_sigma", "alpha", "logalpha"),
                 levels = c("sigma", "logit_sigma", "alpha", "logalpha")),
  x     = c(sigma_hat, x_hat, alpha_hat, y_hat)
)

labels <- as_labeller(c(
  sigma        = expression(sigma),
  logit_sigma  = expression(logit(sigma)),
  alpha        = expression(alpha),
  logalpha     = expression(log(alpha))
))

p <- ggplot(df_all, aes(x = x, y = loglik)) +
  geom_line(size = 0.8) +
  geom_vline(data = vlines_df, aes(xintercept = x), linetype = "dashed",
             color = "cornflowerblue") +
  facet_wrap(~ coord, scales = "free_x", labeller = labels) +
  labs(
    x     = NULL,
    y     = "log EPPF",
    title = "PY(5,0.7) log EPPF and EB estimates"
  ) +
  theme_minimal()

print(p)

est_pics_path <- here("res", "pics", "estimation")
ggsave(filename = paste0(est_pics_path, "/EB_logeppf_1000_PY_07.png"),
       plot     = p,
       dpi      = 300,
       width = 6,
       height = 5,
       bg       = "white")
