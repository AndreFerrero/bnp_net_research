library(here)
library(bayesplot)
library(ggplot2)
library(here)
library(tidyverse)
library(posterior)

poll_dir <- here("pollinators/mpl025")

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

options(
  mc.cores = parallel::detectCores(),
  rstan.auto_write = TRUE
)

stan_folder <- here("code", "stan")
ppc_path <- here(stan_folder, "basic_ppc.stan")
ppc_mod <- stan_model(file = ppc_path)

data <- list(
  K_A = nrow(plant),
  K_B = nrow(poll),
  n_A = plant$counts,
  n_B = poll$counts,
  prior_alpha_A = c(3, 0.4),
  prior_alpha_B = c(3, 0.4),
  prior_sigma_A = c(1.5, 3),
  prior_sigma_B = c(1.5, 3),
  e_obs = sum(edges$weight)
)

# Fit model
fit <- sampling(
  object = ppc_mod,
  data   = data,
  chains = 4,
  iter   = 5000,
  warmup = 1000,
  seed   = 42,
)

save(fit, file = here(poll_dir, "weak_fit.Rdata"))

mcmc_acf(fit)
mcmc_trace(fit)
mcmc_dens(fit)

summary(fit)$summary

draws_df <- as_draws_df(fit)
d_ppc <- draws_df$density_ppc

# Calculate p-value
p_value <- mean(d_ppc >= d_obs)
p_value_label <- paste("Posterior p-value =", format(p_value, digits = 3))

# Calculate line position robustly
prelim_plot <- ggplot(data.frame(density = d_ppc), aes(x = density)) +
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
  labs(
    x = "Density",
    y = NULL
  ) +
  consistent_theme +
  theme(
    axis.ticks.y = element_blank(), # removes y-axis tick marks
    axis.text.y = element_blank(), # removes y-axis labels
    axis.line.y = element_blank(), # removes the axis line if present
    panel.grid.major.y = element_blank(), # removes major horizontal grid lines
    panel.grid.minor.y = element_blank() # removes minor horizontal grid lines
  )

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
