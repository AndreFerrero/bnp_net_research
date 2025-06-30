library(rstan)
library(here)
library(posterior)
library(bayesplot)
library(ggplot2)
library(tidyverse)
library(patchwork)

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


# --- The Updated Function ---
bayes_plots <- function(
    fit_object, color_set, net, pics_folder,
    alpha_true, sigma_true, save = FALSE) {
  # --- 1. Setup: Define colors and a consistent theme for ALL plots ---
  color_scheme_set(color_set)
  current_colors <- color_scheme_get()
  main_color <- current_colors$light
  highlight_color <- current_colors$dark

  param_labels <- c(
    alpha_A = "alpha[A]", alpha_B = "alpha[B]",
    sigma_A = "sigma[A]", sigma_B = "sigma[B]"
  )

  true_values_df <- data.frame(
    parameter = names(param_labels),
    true_value = c(alpha_true, sigma_true)
  )

  draws <- as_draws_df(fit_object)

  # Define a single theme object to apply to all plots for consistency
  consistent_theme <- theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(size = 10),
      plot.margin = margin(t = 5, r = 10, b = 5, l = 5),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "plain", size = 11)
    )

  # --- 2. Diagnostic & Posterior Plots with Harmonized Style ---

  # ACF plot
  acf_plot <- mcmc_acf_bar(fit_object,
    pars = names(param_labels),
    facet_args = list(labeller = ggplot2::labeller(
      .default = label_parsed,
      .cols = param_labels
    ))
  ) +
    consistent_theme +
    hline_at(0.2, linetype = 2, size = 0.15, color = "gray") +
    scale_y_continuous(breaks = c(0.2))

  # Trace plot
  trace_plot <- mcmc_trace(fit_object,
    pars = names(param_labels),
    facet_args = list(
      nrow = 2,
      labeller = ggplot2::labeller(
        .default = label_parsed,
        .cols = param_labels
      )
    )
  ) +
    consistent_theme

  # Posterior densities
  dens_plot <- mcmc_dens_overlay(fit_object,
    pars = names(param_labels),
    facet_args = list(
      nrow = 2,
      labeller = ggplot2::labeller(
        .default = label_parsed,
        .cols = param_labels
      )
    )
  ) +
    consistent_theme +
    scale_y_continuous(labels = NULL)

  make_param_plot <- function(param_name, true_val, color, parsed_title) {
    # Extract draws
    df <- draws %>%
      select({{ param_name }}) %>%
      rename(value = {{ param_name }})

    # Compute CI
    ci <- quantile(df$value, probs = c(0.025, 0.975))
    ci_low <- ci[1]
    ci_high <- ci[2]

    # Expand CI window by ±30%
    ci_width <- ci_high - ci_low
    xmin <- max(0, ci_low - 0.3 * ci_width)
    xmax <- ci_high + 0.3 * ci_width

    # Density
    dens <- density(df$value)
    dens_df <- tibble(x = dens$x, y = dens$y) %>%
      filter(x >= xmin & x <= xmax) %>%
      mutate(in_ci = x >= ci_low & x <= ci_high)

    # Label position (centered along x, low y)
    ci_label_x <- (ci_low + ci_high) / 2
    ci_label_y <- 0 # baseline near x axis ticks

    # Plot
    ggplot() +
      # CI shaded area
      geom_area(
        data = filter(dens_df, in_ci),
        aes(x = x, y = y),
        fill = main_color,
        alpha = 0.7
      ) +
      # Density line
      geom_line(
        data = dens_df,
        aes(x = x, y = y),
        color = color,
        linewidth = 0.8
      ) +
      # True value line (starts from y=0)
      annotate(
        "segment",
        x = true_val, xend = true_val,
        y = 0, yend = max(dens_df$y),
        color = "black",
        linewidth = 0.8,
        linetype = "solid"
      ) +
      # "95% CI" text near x-axis ticks
      annotate(
        "text",
        x = ci_label_x,
        y = 0,
        label = "95% CI",
        vjust = 2,
        size = 3.5,
        color = highlight_color
      ) +
      coord_cartesian(clip = "off") +
      # CI ticks with 2 decimal labels
      scale_x_continuous(
        breaks = c(round(ci_low, 3), round(ci_high, 3)),
        limits = c(xmin, xmax),
        labels = scales::label_number(accuracy = 0.01)
      ) +
      scale_y_continuous(
        breaks = NULL
      ) +
      ggtitle(parsed_title) +
      consistent_theme +
      labs(x = NULL, y = NULL) +
      theme(
        plot.title = element_text(hjust = 0.5)
      )
  }

  plot_alpha_A <- make_param_plot("alpha_A", alpha_true[1], highlight_color, expression(alpha[A]))
  plot_alpha_B <- make_param_plot("alpha_B", alpha_true[2], highlight_color, expression(alpha[B]))
  plot_sigma_A <- make_param_plot("sigma_A", sigma_true[1], highlight_color, expression(sigma[A]))
  plot_sigma_B <- make_param_plot("sigma_B", sigma_true[2], highlight_color, expression(sigma[B]))


  (estimates_plot <- (plot_alpha_A | plot_alpha_B) /
    (plot_sigma_A | plot_sigma_B))

  # --- 3. PPC Plot with p-value annotation ---
  d_obs <- nrow(unique(net$edges)) / (net$xA$K * net$xB$K)
  unif_d_ppc <- draws$density_ppc

  # Calculate p-value
  p_value <- mean(unif_d_ppc >= d_obs)
  p_value_label <- paste("Posterior p-value =", format(p_value, digits = 3))

  # Calculate line position robustly
  prelim_plot <- ggplot(data.frame(density = unif_d_ppc), aes(x = density)) +
    geom_histogram(bins = 30)
  plot_data <- ggplot_build(prelim_plot)
  bin_data <- plot_data$data[[1]]

  # Check if d_obs is within the range of simulations to avoid errors
  target_bin <- bin_data %>% filter(d_obs >= xmin & d_obs < xmax)
  if (nrow(target_bin) > 0) {
    line_position <- target_bin$xmin
  } else {
    line_position <- d_obs
  }

  # --- THE DYNAMIC FIX: Prepare data for ggrepel ---
  # Isolate the data for only the bars in the tail (the highlighted ones)
  tail_bars <- bin_data %>%
    filter(xmin >= line_position)

  # Create a data frame for ggrepel. We will add the label to only ONE row.
  # The other rows will act as invisible "repulsion points".
  if (nrow(tail_bars) > 0) {
    # Find the tallest bar in the tail to anchor our label
    label_anchor_index <- which.max(tail_bars$count) + 1

    # Create a new column for the label text
    tail_bars$label_text <- ""
    tail_bars$label_text[label_anchor_index] <- p_value_label
  }
  # --- End of Fix ---


  ppc_plot <- ggplot(data.frame(density = unif_d_ppc), aes(x = density)) +
    geom_histogram(
      aes(fill = after_stat(xmin) >= line_position),
      bins = 30, color = "white", linewidth = 0.2
    ) +
    scale_fill_manual(values = c("FALSE" = main_color, "TRUE" = highlight_color)) +
    annotate(
      "segment",
      x = line_position, xend = line_position,
      y = 0, yend = Inf,
      color = "black",
      size = 1.6
    ) +
    annotate("text",
      x = line_position, y = Inf, label = paste("Observed:", round(d_obs, 3)),
      vjust = 2, color = "black", hjust = -0.1, fontface = "plain",
      size = 4
    ) +
    labs(
      x = "Density",
      y = "Frequency"
    ) +
    consistent_theme

  # Add the ggrepel layer ONLY if there are tail bars to label
  if (nrow(tail_bars) > 0) {
    ppc_plot <- ppc_plot +
      # Use geom_text_repel instead of annotate
      ggrepel::geom_text_repel(
        data = tail_bars,
        aes(x = x, y = count, label = label_text),
        color = highlight_color,
        size = 4,
        fontface = "plain",
        # --- Repulsion control ---
        box.padding = 2, # How far the label box should be from points/other labels
        point.padding = 0.3, # How far the label box should be from its own anchor point
        min.segment.length = 5
      )
  }

  # --- 4. Saving logic (unchanged) ---
  if (save) {
    ggsave(filename = here(pics_folder, "unif_ppc_acf.pdf"), plot = acf_plot, device = "pdf", width = 4, height = 3)
    ggsave(filename = here(pics_folder, "unif_ppc_trace.pdf"), plot = trace_plot, device = "pdf", width = 4, height = 3)
    ggsave(filename = here(pics_folder, "unif_ppc_post.pdf"), plot = dens_plot, width = 5, height = 3)
    ggsave(filename = here(pics_folder, "unif_ppc_est.pdf"), plot = estimates_plot, width = 5, height = 3)
    ggsave(filename = here(pics_folder, "unif_ppc_dens.pdf"), plot = ppc_plot, width = 5, height = 3)
  }
}

bayes_plots(unif_0_02_ppc_fit, "red", net_0_02,
  unif_0_02_pics_folder, alpha_true, sigma_true_0_02,
  save = TRUE
)
bayes_plots(unif_0_07_ppc_fit, "blue", net_0_07,
  unif_0_07_pics_folder, alpha_true, sigma_true_0_07,
  save = TRUE
)
