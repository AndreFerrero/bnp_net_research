# --- Load Libraries ---
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(here)

# --- Load and Prepare Data (same as before) -------------------------------
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

graph <- graph_from_data_frame(full_edges, vertices = nodes, directed = FALSE)

# --- Colours & node degree (same logic) -----------------------------------
n_plants <- length(plants)
plant_palette <- if (n_plants <= 12) {
  brewer.pal(n_plants, "Paired")
} else {
  colorRampPalette(brewer.pal(12, "Paired"))(n_plants)
}

plant_colors <- setNames(plant_palette, plants)
plant_colors[4] <- "#ED4F50"
plant_colors[9] <- "#2A7FB7"
plant_colors[12] <- "gold2"

V(graph)$color <- ifelse(V(graph)$name %in% names(plant_colors),
  plant_colors[V(graph)$name],
  "gray50"
)
E(graph)$color <- plant_colors[full_edges$plant]
V(graph)$degree <- degree(graph, mode = "all")

# --- Build coordinates ----------------------------------------------------
# layout_as_bipartite puts plants at y = 0 and pollinators at y = 1.
# We keep that, then flip later with coord_flip().
# Get initial layout
layout <- layout_as_bipartite(graph)

# Identify node groups
plant_ids <- which(V(graph)$type == FALSE)
pollinator_ids <- which(V(graph)$type == TRUE)

# Get range of original x positions for plants and pollinators
plant_x_range <- range(layout[plant_ids, 1])
pollinator_x_range <- range(layout[pollinator_ids, 1])

# Assign evenly spaced x-positions within original range
layout[plant_ids, 1] <- seq(from = pollinator_x_range[1], to = pollinator_x_range[2], length.out = length(plant_ids))
layout[pollinator_ids, 1] <- seq(from = pollinator_x_range[1], to = pollinator_x_range[2], length.out = length(pollinator_ids))

# Wrap as data frame
coords <- as.data.frame(layout)
colnames(coords) <- c("x", "y")
coords$name <- V(graph)$name
coords$type <- V(graph)$type
coords$color <- V(graph)$color
coords$degree <- V(graph)$degree

# Edges: join source & target coordinates + edge colour --------------------
edges_df <- full_edges %>%
  left_join(coords, by = c("plant" = "name")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(coords, by = c("pollinator" = "name")) %>%
  rename(x2 = x, y2 = y) %>%
  mutate(color = plant_colors[plant])

# --- Plot -----------------------------------------------------------------
library(ggplot2)
graph_plot <- ggplot() +
  geom_segment(
    data = edges_df,
    aes(
      x = x1, y = y1, xend = x2, yend = y2,
      colour = I(color)
    ),
    alpha = 0.6, linewidth = 0.7
  ) +
  geom_point(
    data = coords,
    aes(x = x, y = y, colour = I(color), size = degree)
  ) +
  # Removed geom_text to hide node names
  scale_size(range = c(2, 6)) +
  theme_void() +
  theme(legend.position = "none")

# Calculate mean x position for each group to center labels horizontally
label_x_plants <- mean(coords$x[coords$type == FALSE])
label_x_pollinators <- mean(coords$x[coords$type == TRUE])

graph_plot <- graph_plot +
  annotate("text",
    x = label_x_plants,
    y = max(coords$y[coords$type == FALSE]) + 0.1, # above plants row
    label = "Plants",
    fontface = "bold",
    size = 5,
    hjust = 0.5 # center horizontally
  ) +
  annotate("text",
    x = label_x_pollinators,
    y = min(coords$y[coords$type == TRUE]) - 0.1, # below pollinators row
    label = "Pollinators",
    fontface = "bold",
    size = 5,
    hjust = 0.5 # center horizontally
  )


# --- Save -----------------------------------------------------------------
ggsave(here(poll_dir, "graph_bipartite_even_h.pdf"),
  graph_plot,
  width = 7, height = 5
)
