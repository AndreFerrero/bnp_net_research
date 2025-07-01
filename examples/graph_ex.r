library(igraph)

set.seed(123)
# --- Setup (same as before) ---
edges <- data.frame(
  from = c("A", "A", "B", "C", "C", "D", "A"),
  to = c("B", "C", "C", "D", "E", "E", "B"),
  label = c("1", "2", "3", "4", "5", "6", "7")
)
g <- graph_from_data_frame(edges, directed = FALSE)
layout <- layout_with_fr(g)
E(g)$curved <- curve_multiple(g)


# --- 1) Save Plain graph (no labels) to PDF ---
pdf("graph_plain.pdf", width = 3, height = 2)
par(mar = c(0, 0, 0, 0)) # Set margins to zero for a clean look
plot(g,
  layout = layout,
  vertex.label = NA,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  edge.label = NA
)
dev.off()

# --- 2) Save Graph with vertex labels to PDF ---
pdf("graph_vertex_labels.pdf", width = 3, height = 2)
par(mar = c(0, 0, 0, 0))
plot(g,
  layout = layout,
  vertex.label = V(g)$name,
  vertex.label.color = "black",
  vertex.label.cex = 1.2,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  edge.label = NA
)
dev.off()

# --- 3) Save Graph with edge labels to PDF ---
pdf("graph_edge_labels.pdf", width = 3, height = 2)
par(mar = c(0, 0, 0, 0))
plot(g,
  layout = layout,
  vertex.label = NA,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  edge.label = E(g)$label,
  edge.label.color = "black",
  edge.label.cex = 1.5
)
dev.off()

pdf("combined_graphs.pdf", width = 4, height = 2)
par(mfrow = c(1, 3), mar = c(0, 0, 2, 0))
plot(g,
  layout = layout,
  vertex.label = NA,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  edge.label = NA,
  main = "(a)"
)

plot(g,
  layout = layout,
  vertex.label = V(g)$name,
  vertex.label.color = "black",
  vertex.label.cex = 1.2,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  edge.label = NA,
  main = "(b)"
)

plot(g,
  layout = layout,
  vertex.label = NA,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  edge.label = E(g)$label,
  edge.label.color = "black",
  edge.label.cex = 1.5,
  main = "(c)"
)

dev.off()
par(mfrow = c(1, 1))
