# Load required package
library(igraph)

# Define the edges (including duplicates)
edges <- c("A", "B", "A", "B", "A", "C")

# Create the graph with multiple edges
g <- graph(edges, directed = FALSE)

# Automatically assign curvature for multiple edges
E(g)$curved <- curve_multiple(g)

# Save layout to keep positions consistent across both plots
layout_fixed <- layout_with_kk(g)

# Open PDF device
pdf("graph_comparison.pdf", width = 12, height = 6) # Wider for side-by-side

# Set up plotting area: 1 row, 2 columns
par(mfrow = c(1, 2))

# ---------- Plot 1: Repeated edges, default node shapes ----------
plot(
  g,
  layout = layout_fixed,
  main = "Repeated Edges (Default Shapes)",
  vertex.size = 30,
  vertex.label.cex = 1.5,
  vertex.color = "lightblue",
  vertex.frame.color = "black"
)

# ---------- Plot 2: Same graph, with different node shapes ----------
# Assign shapes and colors based on group
V(g)$shape <- ifelse(V(g)$name == "A", "circle", "square")
V(g)$color <- ifelse(V(g)$name == "A", "lightblue", "orange")

plot(
  g,
  layout = layout_fixed,
  main = "Node Shapes by Group",
  vertex.size = 30,
  vertex.label.cex = 1.5,
  vertex.frame.color = "black"
)

# Close the PDF device
dev.off()
