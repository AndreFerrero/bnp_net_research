# --- Required Libraries ---
library(rstan)
library(tidyverse)
library(here)
library(cowplot) # For arranging plots
library(ggdist)  # For a nice ggplot theme

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
    legend.position = "top",  # Keep legend on top
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11)
  )

# R-hat Plot (legend on top with parsed labels)
p_rhat <- ggplot(plot_data_agg, aes(x = n_edges_factor, y = max_rhat, color = `Simulation Setup`, group = `Simulation Setup`)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_discrete(labels = function(x) parse(text = x_axis_labels[x])) +  # map factor levels to parsed labels
  scale_y_continuous(name = expression(Maximum~hat(R))) +
  scale_color_viridis_d(
    option = "B",
    end = 0.8,
    labels = parse(text = levels(factor(plot_data_agg$`Simulation Setup`)))
  ) +
  common_theme +
  theme(legend.position = "top")


# Bulk ESS Distribution Plot (no title/subtitle)
p_ess <- ggplot(plot_data_full, aes(x = n_edges_factor, y = Bulk_ESS, fill = `Simulation Setup`)) +
  geom_boxplot(alpha = 0.8) +
  scale_x_discrete(labels = function(x) parse(text = x_axis_labels[x])) +
  scale_y_continuous(name = "Distribution of Bulk ESS") +
  scale_fill_viridis_d(option = "B", end = 0.8, guide = "none") +
  common_theme +
  theme(legend.position = "none")

# Divergent Transitions Plot (no title/subtitle)
p_div <- ggplot(plot_data_agg, aes(x = n_edges_factor, y = total_divergences, color = `Simulation Setup`, group = `Simulation Setup`)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_discrete(name = "Sample Size (Edges)", labels = function(x) parse(text = x_axis_labels[x])) +
  scale_y_continuous(name = "Divergent Transitions", limits = c(0, NA)) +
  scale_color_viridis_d(option = "B", end = 0.8) +
  common_theme +
  theme(axis.title.x = element_text(size = 14),
  legend.position = "none")

# --- Combine and save the final visualization ---

# Arrange the three plots vertically
final_plot <- plot_grid(
  p_rhat,
  p_ess,
  p_div,
  ncol = 1,
  align = 'v',
  rel_heights = c(1.15, 1, 1)
)

# Define output path
output_dir <- here("code", "estimation")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
output_path <- file.path(output_dir, "convergence_visualization.pdf")

# Save the plot
ggsave(output_path, final_plot, width = 9, height = 11, bg = "white")

message(sprintf("\n--- Part 2 complete. Hybrid visualization saved to: ---\n%s", output_path))