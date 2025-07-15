# --- Load Libraries ---
library(igraph)
library(ggraph)
library(tidyverse)
library(RColorBrewer) # Added for better qualitative palettes
# library(viridis) # No longer needed for this solution
library(here)

# --- Load and Prepare Data (No changes here) ---
poll_dir <- here("pollinators/mpl025")
raw_mat <- read.csv(here(poll_dir, "M_PL_025.csv"),
  check.names = FALSE, row.names = 1
)

full_edges <- raw_mat %>%
  rownames_to_column("plant") %>%
  pivot_longer(-plant, names_to = "pollinator", values_to = "weight") %>%
  filter(weight > 0)

plants <- unique(full_edges$plant)
pollinators <- unique(full_edges$pollinator)

nodes <- tibble(
  name = c(plants, pollinators),
  type = c(rep(FALSE, length(plants)), rep(TRUE, length(pollinators)))
)

graph <- graph_from_data_frame(d = full_edges, vertices = nodes, directed = FALSE)

# Color scheme
n_plants <- length(plants)
if (n_plants <= 12) {
  # Use the 'Paired' palette if it has enough colors
  plant_palette <- brewer.pal(n_plants, "Paired")
} else {
  # Otherwise, create a custom palette by interpolating the 'Paired' colors
  plant_palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_plants)
}

# Create a named vector for easy lookup
plant_colors <- setNames(plant_palette, plants)
plant_colors[4] <- "#ED4F50"
plant_colors[9] <- "#2A7FB7"
plant_colors[12] <- "gold2"

# Assign colors to nodes and edges
V(graph)$color <- ifelse(
  V(graph)$name %in% names(plant_colors),
  plant_colors[V(graph)$name],
  "gray50"
)
E(graph)$color <- plant_colors[full_edges$plant]

# Calculate node degree
V(graph)$degree <- degree(graph, mode = "all")


# --- Plot with Improvements ---
graph_plot <- ggraph(graph, layout = "bipartite") +
  geom_edge_link(aes(color = I(color)), alpha = 0.6, width = 0.7) +
  geom_node_point(aes(color = I(color), size = degree)) +

  # Align labels
  geom_node_text(
    aes(label = name),
    # Set alignment: 1 = right align, 0 = left align
    hjust = ifelse(V(graph)$type, 1, 0),
    # Nudge text away from the node center by a fixed amount
    nudge_y = ifelse(V(graph)$type, -0.05, 0.05),
    size = 4
  ) +

  scale_size(range = c(2, 6)) +
  # Expand the horizontal plot limits to make room for the nudged labels
  scale_y_continuous(expand = expansion(mult = c(0.6, 0.65))) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none")

# Save the plot
ggsave(here(poll_dir, "graph_bipartite.pdf"), graph_plot, width = 8, height = 10)
