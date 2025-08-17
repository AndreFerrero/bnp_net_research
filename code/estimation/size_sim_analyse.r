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

# –– 4) MODIFIED Plotting function (with centered titles)
plot_posteriors <- function(df, sigma_lbl) {
  # This function creates four separate plots and combines them using cowplot
  # to allow for fine-grained control over the x-axis of each individual plot.

  # A consistent order for the parameters
  param_order <- c("alpha_A", "alpha_B", "sigma_A", "sigma_B")
  params_to_plot_local <- intersect(param_order, unique(df$parameter))

  plot_list <- list()

  for (param_name in params_to_plot_local) {
    # Filter data for the current parameter
    df_param <- df %>% filter(parameter == param_name)
    
    # Get the unique true value for this parameter
    true_val <- unique(df_param$true_value)
    
    # --- AXIS MODIFICATION START ---
    # 1. Calculate breaks that are guaranteed to include the true value
    breaks <- sort(unique(c(scales::pretty_breaks(n = 3)(df_param$value), true_val)))

    # 2. Use a label function that shows decimals only when necessary
    label_formatter <- scales::label_number_auto()
    # --- AXIS MODIFICATION END ---

    # Create the plot title to mimic the original facet label
    parameter_label <- parse(text = case_when(
        param_name == "alpha_A" ~ "alpha[A]",
        param_name == "alpha_B" ~ "alpha[B]",
        param_name == "sigma_A" ~ "sigma[A]",
        param_name == "sigma_B" ~ "sigma[B]",
        TRUE ~ param_name
    ))

    p_single <- ggplot(df_param, aes(x = value, colour = factor(n_edges))) +
        geom_density(size = 0.5, key_glyph = "path") +
        geom_vline(xintercept = true_val, linetype = "dashed", colour = "red") +
        scale_x_continuous(
            breaks = breaks,
            labels = label_formatter
        ) +
        scale_color_viridis_d(
            option = "B",
            end = 0.9,
            direction = -1,
            name = "Sample size (edges)",
            labels = function(x) parse(text = paste0("10^", log10(as.numeric(x))))
        ) +
        theme_minimal(base_size = 14) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none", # Remove legend from individual plots
            plot.title = element_text(hjust = 0.5) # <<< THIS LINE CENTERS THE TITLE
        ) +
        labs(x = NULL, y = NULL, title = parameter_label)
    
    plot_list[[param_name]] <- p_single
  }
  
  # Arrange the plots into a 2x2 grid.
  p_grid <- plot_grid(plotlist = plot_list, ncol = 2)
  
  return(p_grid)
}

# –– 5) Save one plot per sigma (with minor simplification)
for (lbl in sigma_labels) {
  df_sub <- posterior_df %>% filter(sigma_label == lbl)
  
  # The new function returns a combined plot with no legend,
  # matching the desired output of the original script.
  p <- plot_posteriors(df_sub, lbl)
  
  out_file <- here("code", "estimation", sprintf("subnet_posterior_plot_sigma_%s.pdf", lbl))
  ggsave(out_file, p, width = 10, height = 6, bg = "white")
  message("Saved: ", out_file)
}