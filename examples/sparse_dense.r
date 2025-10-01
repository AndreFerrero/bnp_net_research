# Load igraph
set.seed(26)
library(igraph)

# Set seed for reproducibility
# --- Create Sparse Graph ---
# 10 nodes, 10 edges (low density)
sparse_g <- erdos.renyi.game(n = 9, p.or.m = 10, type = "gnm", directed = FALSE)

# --- Create Dense Graph ---
# 10 nodes, 35 edges (high density, close to complete)
dense_g <- erdos.renyi.game(n = 9, p.or.m = 35, type = "gnm", directed = FALSE)

# --- Common Layout (same for both) ---
layout <- layout_with_fr(sparse_g) # Same layout for fair comparison

# --- Plot Sparse Graph ---
pdf("examples/sparse_graph_black.pdf", width = 5, height = 5)
par(mar = c(0, 0, 0, 0))
plot(sparse_g,
  layout = layout,
  vertex.color = "darkgray",
  vertex.label = NA,
  vertex.size = 20,
  edge.color = "black",
  edge.width = 1
)
dev.off()

# --- Plot Dense Graph ---
pdf("examples/dense_graph_black.pdf", width = 5, height = 5)
par(mar = c(0, 0, 0, 0))
plot(dense_g,
  layout = layout,
  vertex.color = "darkgray",
  vertex.label = NA,
  vertex.size = 20,
  edge.color = "black",
  edge.width = 1
)
dev.off()

pdf("sparse_dense_graphs.pdf", width = 4, height = 2.5)
par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))

plot(dense_g,
  layout = layout,
  vertex.color = "gray90",
  vertex.label = NA,
  vertex.size = 20,
  edge.color = "black",
  edge.width = 1,
  main = "(a)"
)

plot(sparse_g,
  layout = layout,
  vertex.color = "gray90",
  vertex.label = NA,
  vertex.size = 20,
  edge.color = "black",
  edge.width = 1,
  main = "(b)"
)

dev.off()
par(mfrow = c(1, 1))
