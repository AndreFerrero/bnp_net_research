library(tidyverse)
library(here)
library(cowplot) # For arranging plots
library(ggdist)  # For better theme

# --- 1) Load and Prepare the Convergence Data ---

# Define the path to the summary file
summary_csv_path <- here("code", "estimation", "convergence_reports", "convergence_summary.csv")

if (!file.exists(summary_csv_path)) {
  stop("The convergence summary file was not found. Please run the convergence checking script first.\nMissing file: ", summary_csv_path)
}

message("Loading convergence data from: ", summary_csv_path)
convergence_df <- read_csv(summary_csv_path, show_col_types = FALSE)

# --- 2) Prepare Data for Plotting ---

# Add descriptive columns and a factor for the x-axis
plot_data_full <- convergence_df %>%
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

message("Data prepared for visualization.")

# --- 3) Create the Three Diagnostic Plots ---

# Common theme for all plots
common_theme <- theme_minimal_hgrid(14) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(), # Remove x-axis title from individual plots
    axis.text.x = element_text(size = 11)
  )

# Custom scientific labels for the x-axis
x_axis_labels <- c(
  "100" = expression(10^2),
  "1000" = expression(10^3),
  "10000" = expression(10^4),
  "100000" = expression(10^5)
)

# A) R-hat Plot (Maximum Value Line Plot)
p_rhat <- ggplot(plot_data_agg, aes(x = n_edges_factor, y = max_rhat, color = `Simulation Setup`, group = `Simulation Setup`)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_discrete(labels = x_axis_labels) +
  scale_y_continuous(name = expression(Maximum~hat(R))) +
  scale_color_viridis_d(
    option = "B",
    end = 0.8,
    labels = parse(text = levels(factor(plot_data_agg$`Simulation Setup`)))
  ) +
  labs(
    title = "Convergence Diagnostics by Sample Size",
    subtitle = "Worst-case R-hat across all parameters"
  ) +
  common_theme

# B) Bulk ESS Distribution Plot (Boxplots, Linear Scale)
p_ess <- ggplot(plot_data_full, aes(x = n_edges_factor, y = Bulk_ESS, fill = `Simulation Setup`)) +
  geom_boxplot(alpha = 0.8) +
  scale_x_discrete(labels = x_axis_labels) +
  scale_y_continuous(name = "Distribution of Bulk ESS (Linear Scale)") +
  scale_fill_viridis_d(option = "B", end = 0.8, guide = "none") + # No legend
  labs(subtitle = "Effective Sample Size for all parameters") +
  common_theme

# C) Divergent Transitions Plot (Line/Point)
p_div <- ggplot(plot_data_agg, aes(x = n_edges_factor, y = total_divergences, color = `Simulation Setup`, group = `Simulation Setup`)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_discrete(
    name = "Sample Size (Edges)",
    labels = x_axis_labels
  ) +
  scale_y_continuous(
    name = "Divergent Transitions",
    limits = c(0, NA) # Ensure y-axis starts at 0
  ) +
  scale_color_viridis_d(option = "B", end = 0.8, guide = "none") + # No legend
  labs(subtitle = "Count of divergent transitions post-warmup") +
  common_theme +
  theme(axis.title.x = element_text(size = 14)) # Add back x-axis title for the bottom plot

# --- 4) Combine Plots into a Single Visualization ---

# Arrange the three plots vertically
final_plot <- plot_grid(
  p_rhat,
  p_ess,
  p_div,
  ncol = 1,
  align = 'v', # Align them vertically
  rel_heights = c(1.15, 1, 1) # Give slightly more space to the top plot with the legend
)

# --- 5) Save the Final Plot ---
output_dir <- here("code", "estimation")
output_path <- file.path(output_dir, "convergence_visualization.pdf")

ggsave(output_path, final_plot, width = 9, height = 11, bg = "white")

message("\nHybrid visualization complete! Plot saved to:\n", output_path)