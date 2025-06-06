library(here)
library(rstan)
library(dplyr)

est_f <- here("code", "estimation")

true_vals <- data.frame(
  Parameter = c("alpha_A", "alpha_B", "sigma_A", "sigma_B"),
  True_1 = c(5, 5, 0, 0.7),
  True_2 = c(5, 5, 0, 0.2)
)

load(here(est_f, "0_07_unif_ppc_fit.Rdata"))
fit1 <- unif_ppc_fit
rm(unif_ppc_fit)

load(here(est_f, "0_02_unif_ppc_fit.Rdata"))
fit2 <- unif_ppc_fit
rm(unif_ppc_fit)

# Function to extract summary for a fit
extract_summary <- function(fit, params) {
  summ <- summary(fit, pars = params)$summary
  data.frame(
    Mean = round(summ[, "mean"], 2),
    CI = paste0("[", round(summ[, "2.5%"], 2), ", ", round(summ[, "97.5%"], 2), "]")
  )
}

# Parameters to extract
params <- c("alpha_A", "alpha_B", "sigma_A", "sigma_B")

# Extract summaries
summ1 <- extract_summary(fit1, params)
summ2 <- extract_summary(fit2, params)

# Combine into one data frame
latex_df <- cbind(true_vals, summ1, summ2)
colnames(latex_df) <- c(
  "Parameter", "True_1", "True_2",
  "Mean_1", "CI_1", "Mean_2", "CI_2"
)

# Replace plain names with LaTeX math symbols
latex_df$Parameter <- recode(latex_df$Parameter,
  "alpha_A" = "\\alpha_A",
  "alpha_B" = "\\alpha_B",
  "sigma_A" = "\\sigma_A",
  "sigma_B" = "\\sigma_B"
)

# Print LaTeX table
cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{Posterior summaries}\n")
cat("\\label{tab:posterior_summary}\n")
cat("\\begin{tabular}{lcccccc}\n")
cat("\\toprule\n")
cat("Parameter & True (1) & Mean (1) & 95\\% CI (1) & True (2) & Mean (2) & 95\\% CI (2) \\\\\n")
cat("\\midrule\n")

for (i in 1:nrow(latex_df)) {
  cat(paste0(
    "$", latex_df$Parameter[i], "$ & ",
    latex_df$True_1[i], " & ", latex_df$Mean_1[i], " & ", latex_df$CI_1[i], " & ",
    latex_df$True_2[i], " & ", latex_df$Mean_2[i], " & ", latex_df$CI_2[i], " \\\\\n"
  ))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")

# Extract diagnostic summary
extract_diagnostics <- function(fit, params) {
  summ <- summary(fit, pars = params)$summary
  data.frame(
    Rhat = round(summ[, "Rhat"], 3),
    n_eff = round(summ[, "n_eff"], 0)
  )
}

# Extract diagnostics
diag1 <- extract_diagnostics(fit1, params)
diag2 <- extract_diagnostics(fit2, params)

# Get number of divergences
get_divergences <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
}

div1 <- get_divergences(fit1)
div2 <- get_divergences(fit2)

# Combine into one data frame
diag_df <- data.frame(
  Parameter = params,
  Rhat_1 = diag1$Rhat,
  Neff_1 = diag1$n_eff,
  Rhat_2 = diag2$Rhat,
  Neff_2 = diag2$n_eff
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
cat("\\caption{Convergence diagnostics for two simulation setups, including $\\widehat{R}$, effective sample size $n_{\\mathrm{eff}}$, and number of divergent transitions.}\n")
cat("\\label{tab:diagnostics}\n")
cat("\\begin{tabular}{lcccc}\n")
cat("\\toprule\n")
cat("Parameter & $\\widehat{R}$ (1) & $n_{\\mathrm{eff}}$ (1) & $\\widehat{R}$ (2) & $n_{\\mathrm{eff}}$ (2) \\\\\n")
cat("\\midrule\n")

for (i in 1:nrow(diag_df)) {
  cat(paste0(
    "$", diag_df$Parameter[i], "$ & ",
    diag_df$Rhat_1[i], " & ", diag_df$Neff_1[i], " & ",
    diag_df$Rhat_2[i], " & ", diag_df$Neff_2[i], " \\\\\n"
  ))
}

# Add divergences row
cat("\\midrule\n")
cat(paste0("Divergences & \\multicolumn{2}{c}{", div1, "} & \\multicolumn{2}{c}{", div2, "} \\\\\n"))

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
