library(igraph)

set.seed(123)

# --- Create the original graph ---
edges <- data.frame(
  from = c("A", "A", "A"),
  to   = c("B", "C", "D")
)
g1 <- graph_from_data_frame(edges, directed = FALSE)

# --- Generate a permutation of vertex labels ---
original_labels <- V(g1)$name
permuted_labels <- c("C", "B", "A", "D")
name_map <- setNames(permuted_labels, original_labels)

# --- Create permuted graph g2 ---
g2 <- g1
V(g2)$name <- name_map[V(g1)$name]

# --- Use the same layout for both ---
layout <- layout_with_fr(g1)

# --- Determine which labels were permuted ---
# This compares positions, not just names
changed_labels <- V(g1)$name != V(g2)$name
label_colors_g2 <- ifelse(changed_labels, "red", "black")

# --- Plot both graphs ---
pdf("examples/node_exch.pdf", width = 4, height = 2.5)
par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))

# Original graph
plot(g1,
  layout = layout,
  vertex.label = V(g1)$name,
  vertex.label.color = "black",
  vertex.color = "darkgrey",
  vertex.size = 30,
  edge.color = "black",
  main = "(a)"
)

# Permuted graph with highlighted labels
plot(g2,
  layout = layout,
  vertex.label = V(g2)$name,
  vertex.label.color = label_colors_g2,
  vertex.color = "darkgrey",
  vertex.size = 30,
  edge.color = "black",
  main = "(b)"
)

dev.off()
