# Load igraph
library(igraph)

# Define edge list with weights
edges <- data.frame(
  from = c("A_1", "A_1", "A_2", "A_3", "A_3"),
  to = c("B_1", "B_2", "B_3", "B_1", "B_3"),
  weight = c(1, 3, 2, 5, 4) # example weights
)

# Create graph
g <- graph_from_data_frame(edges, directed = FALSE)

# Identify bipartite node types
V(g)$type <- grepl("^A_", V(g)$name)

# Assign shapes and colors
V(g)$shape <- ifelse(V(g)$type, "circle", "square")
V(g)$color <- ifelse(V(g)$type, "lightblue", "orange")

# Create subscript labels using plotmath
labels <- c(
  expression(A[1]), expression(A[2]), expression(A[3]),
  expression(B[1]), expression(B[2]), expression(B[3])
)

# Compute bipartite layout
layout <- layout_as_bipartite(g)
layout[, 2] <- -layout[, 2] # Flip to have A nodes on top

# Normalize edge weights and assign red scale colors
weights <- E(g)$weight
norm_weights <- (weights - min(weights)) / (max(weights) - min(weights))
edge_colors <- rgb(1, 1 - norm_weights, 1 - norm_weights)

# Plot with subscripts
pdf("bip_plot.pdf", width = 7, height = 4)
plot(g,
  layout = layout,
  vertex.label = labels,
  vertex.label.color = "black",
  vertex.label.cex = 1.4,
  vertex.label.family = "sans",
  vertex.size = 50,
  vertex.shape = V(g)$shape,
  vertex.color = V(g)$color,
  edge.width = weights,
  edge.color = edge_colors,
)

dev.off()
