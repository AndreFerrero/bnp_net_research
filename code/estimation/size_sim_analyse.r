library(rstan)
library(tidyverse)
library(here)
library(bayesplot)
library(cowplot)

# Edge sizes and sigma labels
edge_sizes <- c(1e2, 1e3, 1e4, 1e5)
sigma_labels <- c("0_02", "0_07")

# Parameters to extract
params_to_plot <- c("alpha_A", "alpha_B", "sigma_A", "sigma_B")

# –– 1) True‐value lookup table
true_lookup <- tribble(
  ~sigma_label, ~parameter,   ~true_value,
  "0_02",       "alpha_A",    5,
  "0_02",       "alpha_B",    5,
  "0_02",       "sigma_A",    0.0,
  "0_02",       "sigma_B",    0.2,
  "0_07",       "alpha_A",    5,
  "0_07",       "alpha_B",    5,
  "0_07",       "sigma_A",    0.0,
  "0_07",       "sigma_B",    0.7
)

# –– 2) Read in all fits and stack posteriors
posterior_data <- list()

for (sigma_label in sigma_labels) {
  for (n in edge_sizes) {
    fname <- sprintf("fit_subnet_n_%g_sigma_%s.Rdata", n, sigma_label)
    fpath <- here("code", "estimation", fname)
    if (!file.exists(fpath)) next

    load(fpath) # loads object 'fit'

    post_df <- as.data.frame(
      rstan::extract(fit, pars = params_to_plot)
    ) %>%
      pivot_longer(everything(),
        names_to = "parameter",
        values_to = "value"
      ) %>%
      mutate(
        n_edges     = n,
        sigma_label = sigma_label
      )

    posterior_data[[length(posterior_data) + 1]] <- post_df
  }
}

posterior_df <- bind_rows(posterior_data)

# –– 3) Join true values
posterior_df <- posterior_df %>%
  left_join(true_lookup, by = c("sigma_label", "parameter"))

# –– 4) Plotting function (no title)
plot_posteriors <- function(df, sigma_lbl) {
  df <- df %>%
    mutate(
      parameter_label = case_when(
        parameter == "alpha_A" ~ "alpha[A]",
        parameter == "alpha_B" ~ "alpha[B]",
        parameter == "sigma_A" ~ "sigma[A]",
        parameter == "sigma_B" ~ "sigma[B]",
        TRUE ~ parameter
      )
    )

  ggplot(df, aes(x = value, colour = factor(n_edges))) +
    geom_density(size = 0.5, key_glyph = "path") +  # <- force legend as lines
    geom_vline(aes(xintercept = true_value),
      linetype = "dashed", colour = "red"
    ) +
    facet_wrap(
      ~parameter_label,
      scales = "free",
      ncol = 2,
      labeller = label_parsed
    ) +
    scale_color_viridis_d(
      option = "B",
      end = 0.9,
      direction = -1,
      name = "Sample size (edges)",
      labels = function(x) parse(text = paste0("10^", log10(as.numeric(x)))), # scientific labels
      guide = guide_legend(
        override.aes = list(linetype = 1, shape = NA, fill = NA)  # line only
      )
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(x = NULL, y = NULL, title = NULL)
}

# –– 5) Save one plot per sigma
for (lbl in sigma_labels) {
  df_sub <- posterior_df %>% filter(sigma_label == lbl)
  p <- plot_posteriors(df_sub, lbl)

  out_file <- here("code", "estimation", sprintf("subnet_posterior_plot_sigma_%s.pdf", lbl))
  ggsave(out_file, p, width = 10, height = 6, bg = "white")
  message("Saved: ", out_file)
}

library(patchwork)

# –– 6) Create the two plots separately
plots <- list()
for (lbl in sigma_labels) {
  df_sub <- posterior_df %>% filter(sigma_label == lbl)
  plots[[lbl]] <- plot_posteriors(df_sub, lbl)
}

# –– 7) Combine side by side with shared legend
combined_plot <- (plots[[1]] + plots[[2]]) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box.just = "center"
  )

# –– 11) Save the combined PDF
out_file <- here("code", "estimation", "subnet_sim_both.pdf")
ggsave(out_file, combined_plot, width = 7, height = 4, bg = "white")
