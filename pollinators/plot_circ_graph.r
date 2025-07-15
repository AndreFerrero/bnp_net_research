library(rstan)
library(bayesplot)
library(ggplot2)
library(here)
library(igraph)
library(ggraph)
library(tidyverse)

# Set up directories
poll_dir <- here("pollinators/mpl025")

# Load full data
raw_mat <- read.csv(here(poll_dir, "M_PL_025.csv"),
  check.names = FALSE, row.names = 1
)

# Convert matrix to edge list
full_edges <- raw_mat %>%
  rownames_to_column("plant") %>%
  pivot_longer(-plant, names_to = "pollinator", values_to = "weight") %>%
  filter(weight > 0)

# Create node list
plants <- unique(full_edges$plant)
pollinators <- unique(full_edges$pollinator)

nodes <- tibble(
  name = c(plants, pollinators),
  type = c(rep(TRUE, length(plants)), rep(FALSE, length(pollinators))) # TRUE = plant, FALSE = pollinator
)

# Create graph
graph <- graph_from_data_frame(d = full_edges, vertices = nodes, directed = FALSE)

# Add degree
V(graph)$degree <- degree(graph, mode = "all")

# Create layout: concentric circles
vertex_df <- tibble(
  name = V(graph)$name,
  type = V(graph)$type,
  degree = V(graph)$degree,
  id = 1:vcount(graph)
)

# Split nodes
plants_df <- vertex_df %>%
  filter(type) %>%
  arrange(desc(degree))
pollinators_df <- vertex_df %>%
  filter(!type) %>%
  arrange(desc(degree))

# Radii
r_inner <- 1
r_outer <- 2

# Angles
plants_df <- plants_df %>%
  mutate(
    angle = seq(0, 2 * pi, length.out = n() + 1)[-(n() + 1)],
    x = r_inner * cos(angle),
    y = r_inner * sin(angle)
  )

pollinators_df <- pollinators_df %>%
  mutate(
    angle = seq(0, 2 * pi, length.out = n() + 1)[-(n() + 1)],
    x = r_outer * cos(angle),
    y = r_outer * sin(angle)
  )

# Combine layout by original node order
layout_df <- bind_rows(plants_df, pollinators_df) %>%
  arrange(id) %>%
  select(x, y)

# Plot
graph_plot <- ggraph(graph,
  layout = "manual",
  x = layout_df$x, # <‑‑ x coordinates
  y = layout_df$y
) + # <‑‑ y coordinates
  geom_edge_link(alpha = 0.3, color = "gray50") +
  geom_node_point(aes(color = type, size = degree)) +
  scale_color_manual(
    values = c("black", "forestgreen"),
    labels = c("Pollinator", "Plant")
  ) +
  scale_size(range = c(2, 8)) +
  guides(size = "none") +
  theme_void() +
  theme(legend.position = "top") +
  labs(color = "Node Type")

# Save plot
ggsave(filename = here(poll_dir, "circ_graph.pdf"), graph_plot, width = 7, height = 4)
