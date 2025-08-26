# --- Required Libraries ---
library(rstan)
library(tidyverse)
library(here)
library(cowplot) # For arranging plots
library(ggdist) # For a nice ggplot theme

# ===================================================================
# Part 1: Generate Convergence Summary from .Rdata files
# ===================================================================

message("--- Starting Part 1: Generating convergence summary ---")

# --- Setup for data loading ---
edge_sizes <- c(1e2, 1e3, 1e4, 1e5)
sigma_labels <- c("0_02", "0_07")

# List to store summary data from each model
convergence_summary_list <- list()

# --- Loop through all model fits ---
for (sigma_label in sigma_labels) {
  for (n in edge_sizes) {
    # Define file path
    fname <- sprintf("fit_subnet_n_%g_sigma_%s.Rdata", n, sigma_label)
    fpath <- here("code", "estimation", fname)

    if (!file.exists(fpath)) {
      warning(sprintf("File not found, skipping: %s", fname))
      next
    }

    # Load the 'fit' object from the .Rdata file
    load(fpath)

    # Extract numerical diagnostics using rstan::monitor
    monitor_df <- as.data.frame(rstan::monitor(fit, print = FALSE)) %>%
      rownames_to_column("parameter")

    # Count divergent transitions post-warmup
    sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
    divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

    # Add metadata to the summary
    monitor_df <- monitor_df %>%
      mutate(
        n_edges = n,
        sigma_label = sigma_label,
        divergent_transitions = divergences
      )

    # Store the summary for this model fit in the list
    convergence_summary_list[[length(convergence_summary_list) + 1]] <- monitor_df

    message(sprintf("Processed: %s", fname))
  }
}

# --- Combine all summaries into a single data frame ---
convergence_df <- bind_rows(convergence_summary_list)

message("\n--- Part 1 complete. Convergence data frame created in memory. ---")


# ===================================================================
# Part 2: Create the Convergence Visualization
# ===================================================================

message("\n--- Starting Part 2: Creating the hybrid visualization ---")

# --- Prepare the data for plotting ---

# Add descriptive columns and a factor for the x-axis to the full dataset
plot_data_full <- convergence_df %>%
  filter(parameter %in% c("alpha_A", "alpha_B", "sigma_A", "sigma_B")) %>%
  mutate(
    `Simulation Setup` = case_when(
      sigma_label == "0_02" ~ "sigma[B] == 0.2",
      sigma_label == "0_07" ~ "sigma[B] == 0.7"
    ),
    n_edges_factor = factor(n_edges, levels = c(1e2, 1e3, 1e4, 1e5))
  )

# Create a summarized dataset for the R-hat and Divergence line plots
plot_data_agg <- plot_data_full %>%
  group_by(n_edges_factor, `Simulation Setup`) %>%
  summarise(
    max_rhat = max(Rhat, na.rm = TRUE),
    total_divergences = first(divergent_transitions),
    .groups = "drop"
  )

message("Data prepared for plotting.")

# --- Define custom colors (lighter viridis pair) ---
custom_colors <- c("sigma[B] == 0.2" = "#4D4B4B",
                   "sigma[B] == 0.7" = "#EEAC1D")

# --- Create the three diagnostic plots ---

x_axis_labels <- c(
  "100" = "10^2",
  "1000" = "10^3",
  "10000" = "10^4",
  "1e+05" = "10^5"
)

# Common theme without titles/subtitles
common_theme <- theme_minimal_hgrid(14) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11)
  )

# Max R-hat plot without legend
p_rhat <- ggplot(plot_data_agg, aes(x = n_edges_factor, y = max_rhat,
                                    color = `Simulation Setup`, group = `Simulation Setup`)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_discrete(labels = function(x) parse(text = x_axis_labels[x])) +
  scale_y_continuous(name = expression(Maximum ~ hat(R))) +
  scale_color_manual(values = custom_colors) +
  common_theme +
  theme(legend.position = "none")

# ESS distribution plot without legend
p_ess <- ggplot(plot_data_full, aes(x = n_edges_factor, y = n_eff,
                                    fill = `Simulation Setup`)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  scale_x_discrete(labels = function(x) parse(text = x_axis_labels[x])) +
  scale_y_continuous(name = "Distribution of ESS", limits = c(6000, NA)) +
  scale_fill_manual(values = custom_colors, guide = "none") + # ensures no legend
  common_theme +
  coord_cartesian(ylim = c(6000, NA), clip = "off") +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none"  # removes legend
  )

# Divergent Transitions Plot
p_div <- ggplot(plot_data_agg, aes(x = n_edges_factor, y = total_divergences,
                                   color = `Simulation Setup`, group = `Simulation Setup`)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_discrete(name = "Sample Size (Edges)",
                   labels = function(x) parse(text = x_axis_labels[x])) +
  scale_y_continuous(name = "Divergent Transitions", limits = c(0, NA)) +
  scale_color_manual(values = custom_colors) +
  common_theme +
  theme(
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

# --- Compute ESS ratio relative to total posterior draws ---
# Compute total draws
# Compute total draws
total_draws <- (fit@sim$iter - fit@sim$warmup) * fit@sim$chains

# Define the ESS values you want ticks for
left_ticks <- c(6000, 8000, 10000, 12000)

# Calculate corresponding ratios
right_ticks <- left_ticks / total_draws

# Create plot with secondary axis
p_ess_ratio <- ggplot(plot_data_full, aes(x = n_edges_factor, y = n_eff,
                                                fill = `Simulation Setup`)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
  scale_x_discrete(labels = function(x) parse(text = x_axis_labels[x])) +
  scale_y_continuous(
    name = "ESS",
    limits = c(6000, NA),
    breaks = left_ticks,                           # left axis ticks
    sec.axis = sec_axis(~ . / total_draws,         # map left to right
                        name = "ESS / Total Draws",
                        breaks = right_ticks,     # right axis ticks
                        labels = scales::percent(right_ticks)) # format as %
  ) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  common_theme +
  theme(
    axis.title.y.left = element_text(size = 12),
    axis.title.y.right = element_text(size = 12),
    legend.position = "none"
  )


# --- Combine and save the final visualization ---
library(patchwork)

# Combine the plots horizontally
final_plot <- p_rhat + p_ess +
  plot_layout(ncol = 2, guides = "collect") &   # collect all legends into one
  theme(legend.position = "top")

# Define output path
output_dir <- here("code", "estimation")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_path <- file.path(output_dir, "convergence_visualization.pdf")

# Save the plot
ggsave(output_path, final_plot, width = 11, height = 5, bg = "white")

message(sprintf("\n--- Part 2 complete. Hybrid visualization saved to: ---\n%s", output_path))

# Save Max R-hat plot
ggsave(
  filename = here(output_dir, "max_rhat_plot.pdf"),
  plot = p_rhat,
  width = 6, height = 4
)

# Save ESS plot
ggsave(
  filename = here(output_dir, "ess_plot.pdf"),
  plot = p_ess,
  width = 6, height = 4
)

ggsave(file.path(output_dir, "ess_ratio_plot.pdf"), plot = p_ess_ratio, width = 6, height = 4)
