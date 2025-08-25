library(tidyverse)
library(posterior)
library(bayesplot)
library(ggrepel)
library(ggplot2)
library(here)
library(patchwork)

poll_dir <- here("pollinators/mpl025")

wo_repl_dir <- here(poll_dir, "wo_repl")
plots_dir <- here(wo_repl_dir, "unif_plots")
fits_dir <- here(wo_repl_dir, "unif_fit")

# IMPORT DATA
# Read the CSV file
# (Change the file path to where you've saved your file)
data <- read.csv(here(poll_dir, "M_PL_025.csv"),
  check.names = FALSE, row.names = 1
)

# Convert the adj matrix to an edge list
edges <- data |>
  rownames_to_column(var = "plant") |>
  pivot_longer(
    cols = -plant,
    names_to = "pollinator",
    values_to = "weight"
  ) |>
  filter(weight > 0)

plant <- edges |>
  group_by(plant) |>
  summarise(counts = sum(weight))

poll <- edges |>
  group_by(pollinator) |>
  summarise(counts = sum(weight))

(d_obs <- nrow(edges) / (nrow(plant) * nrow(poll)))


(total_weight <- sum(edges$weight))
(binary_weight <- nrow(edges))

# ESTIMATION

options(
  mc.cores = parallel::detectCores(),
  rstan.auto_write = TRUE
)

stan_folder <- here("code", "stan")
unif_ppc_path <- here(stan_folder, "unif_ppc.stan")
unif_ppc_mod <- stan_model(file = unif_ppc_path)

unif_ppc_data <- list(
  K_A = nrow(plant),
  K_B = nrow(poll),
  n_A = plant$counts,
  n_B = poll$counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  eps = 0.1,
  e_obs = sum(edges$weight)
)

# Fit model
unif_ppc_fit <- sampling(
  object = unif_ppc_mod,
  data = unif_ppc_data,
  chains = 4,
  iter = 20000,
  warmup = 1000,
  seed = 42,
  thin = 10,
  control = list(adapt_delta = 0.999)
)

# save(unif_ppc_fit, file = here(poll_dir, "full_unif_ppc_fit.Rdata"))

# --- LOAD MODEL FIT ---
load(here(fits_dir, "unif_ppc_fit_100pct.Rdata"))
summary(fit)$summary

draws_df <- as_draws_df(fit)
d_ppc <- draws_df$density_ppc

# --- Posterior predictive check ---
p_value <- mean(d_ppc >= d_obs)
p_value_label <- paste("Posterior p-value =", format(p_value, digits = 3))

# Histogram data for PPC
prelim_plot <- ggplot(data.frame(density = d_ppc), aes(x = density)) +
  geom_histogram(bins = 30)
plot_data <- ggplot_build(prelim_plot)
bin_data <- plot_data$data[[1]]

target_bin <- bin_data %>% filter(d_obs >= xmin & d_obs < xmax)
if (nrow(target_bin) > 0) {
  line_position <- target_bin$xmin
} else {
  line_position <- d_obs
}

tail_bars <- bin_data %>%
  filter(xmin >= line_position)
if (nrow(tail_bars) > 0) {
  label_anchor_index <- which.max(tail_bars$count) + 1
  tail_bars$label_text <- ""
  tail_bars$label_text[label_anchor_index] <- p_value_label
}

# --- Theme and colors ---
consistent_theme <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "plain"),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "plain", size = 11)
  )

color_scheme_set("orange")
current_colors <- color_scheme_get()
main_color <- current_colors$light
highlight_color <- current_colors$dark

# --- PPC plot ---
ppc_plot <- ggplot(data.frame(density = d_ppc), aes(x = density)) +
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
    vjust = 2, color = "black", hjust = 1.1, fontface = "plain",
    size = 4
  ) +
  labs(x = "Density", y = NULL) +
  consistent_theme +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

if (nrow(tail_bars) > 0) {
  ppc_plot <- ppc_plot +
    ggrepel::geom_text_repel(
      data = tail_bars,
      aes(x = x, y = count, label = label_text),
      color = highlight_color,
      size = 4,
      fontface = "plain",
      box.padding = 2,
      point.padding = 0.3,
      min.segment.length = 5
    )
}

ggsave(here(plots_dir, "full_ppc.pdf"), ppc_plot, width = 7, height = 4)

# --- Parameter labels ---
param_labels <- c(
  alpha_A = "alpha[A]", alpha_B = "alpha[B]",
  sigma_A = "sigma[A]", sigma_B = "sigma[B]"
)

# --- ACF plot ---
acf_plot <- mcmc_acf(fit,
  pars = names(param_labels),
  facet_args = list(labeller = ggplot2::labeller(
    .default = label_parsed,
    .cols = param_labels
  ))
) +
  consistent_theme +
  hline_at(0.2, linetype = 2, size = 0.15, color = "gray") +
  scale_y_continuous(breaks = c(0.2))

ggsave(here(plots_dir, "full_acf.pdf"), acf_plot, width = 7, height = 4)

# --- Trace plot ---
trace_plot <- mcmc_trace(fit,
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

ggsave(here(plots_dir, "full_trace.pdf"), trace_plot, width = 7, height = 4)

# --- Posterior density plotting function ---
make_posterior_plot <- function(param_name, parsed_label) {
  df <- draws_df %>%
    select({{ param_name }}) %>%
    rename(value = {{ param_name }})

  ci <- quantile(df$value, probs = c(0.025, 0.975))
  ci_low <- ci[1]; ci_high <- ci[2]
  ci_width <- ci_high - ci_low
  xmin <- max(0, ci_low - 0.3 * ci_width)
  xmax <- ci_high + 0.3 * ci_width

  dens <- density(df$value)
  dens_df <- tibble(x = dens$x, y = dens$y) %>%
    filter(x >= xmin & x <= xmax) %>%
    mutate(in_ci = x >= ci_low & x <= ci_high)

  ci_label_x <- (ci_low + ci_high) / 2

  ggplot() +
    geom_area(
      data = filter(dens_df, in_ci),
      aes(x = x, y = y),
      fill = main_color,
      alpha = 0.7
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = y),
      color = highlight_color,
      linewidth = 0.8
    ) +
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
    scale_x_continuous(
      breaks = c(round(ci_low, 3), round(ci_high, 3)),
      limits = c(xmin, xmax),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_y_continuous(breaks = NULL) +
    consistent_theme +
    labs(x = parsed_label, y = NULL) +  # x-axis label now uses parsed text
    theme(
      plot.title = element_blank(),    # remove plot title
      axis.title.x = element_text(face = "plain")  # style x-axis label
    )
}

# --- Create individual posterior plots ---
plot_alpha_A <- make_posterior_plot("alpha_A", expression(alpha[A]))
plot_alpha_B <- make_posterior_plot("alpha_B", expression(alpha[B]))
plot_sigma_A <- make_posterior_plot("sigma_A", expression(sigma[A]))
plot_sigma_B <- make_posterior_plot("sigma_B", expression(sigma[B]))


# --- Combine into single figure ---
estimates_plot <- (plot_alpha_A | plot_alpha_B) / (plot_sigma_A | plot_sigma_B)
ggsave(here(plots_dir, "posterior_estimates.pdf"), estimates_plot, width = 7, height = 4)



