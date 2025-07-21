library(igraph)

set.seed(123)

# Create the graph data
edges <- data.frame(
  from = c("A", "A", "A"),
  to   = c("B", "C", "D")
)
g1 <- graph_from_data_frame(edges, directed = FALSE)

# Create a second graph with permuted labels
original_labels <- V(g1)$name
permuted_labels <- c("C", "B", "A", "D")
name_map <- setNames(permuted_labels, original_labels)

g2 <- g1
V(g2)$name <- name_map[V(g1)$name]

# Use a consistent layout for both graphs
layout <- layout_with_fr(g1)

# Determine which labels changed to color them red in the second plot
changed_labels <- V(g1)$name != V(g2)$name
label_colors_g2 <- ifelse(changed_labels, "red", "black")

# --- Plotting Section with Corrected Margins ---

# Create a directory for the output if it doesn't exist
if (!dir.exists("examples")) {
  dir.create("examples")
}

# Plot original graph with zero margins
pdf("examples/node_exch_g1.pdf", width = 4, height = 2.5)
# Set margins to zero (bottom, left, top, right)
par(mar = c(0, 0, 0, 0))
plot(g1,
  layout = layout,
  vertex.label = V(g1)$name,
  vertex.label.color = "black",
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  main = NA
)
dev.off()

# Plot permuted graph with zero margins
pdf("examples/node_exch_g2.pdf", width = 4, height = 2.5)
# Set margins to zero (bottom, left, top, right)
par(mar = c(0, 0, 0, 0))
plot(g2,
  layout = layout,
  vertex.label = V(g2)$name,
  vertex.label.color = label_colors_g2,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  main = NA
)
dev.off()
