# Load required package
library(igraph)

# Updated edges for unipartite graph (including B–C)
edges1 <- c("A", "B", "A", "B", "A", "C", "B", "C")
g1 <- graph(edges1, directed = FALSE)

# Bipartite graph: relabel B, C to B_1 and B_2
edges2 <- c("A", "B_1", "A", "B_1", "A", "B_2")
g2 <- graph(edges2, directed = FALSE)

# Curved edges
E(g1)$curved <- curve_multiple(g1)
E(g2)$curved <- curve_multiple(g2)

# Define layout (V shape, A at bottom)
layout_v1 <- matrix(c(0, -1, -1, 1, 1, 1), ncol = 2, byrow = TRUE) # A, B, C
layout_v1 <- layout_v1[match(V(g1)$name, c("A", "B", "C")), ]

layout_v2 <- matrix(c(0, -1, -1, 1, 1, 1), ncol = 2, byrow = TRUE) # A, B_1, B_2
layout_v2 <- layout_v2[match(V(g2)$name, c("A", "B_1", "B_2")), ]

# Open PDF
pdf("examples/uni_vs_bip.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))

# Plot 1: Unipartite with B–C edge
plot(
  g1,
  layout = layout_v1,
  main = "Unipartite",
  vertex.size = 50,
  vertex.label.cex = 1.5,
  vertex.color = "lightblue",
  vertex.frame.color = "black"
)

# Plot 2: Bipartite with B[1], B[2] and shapes
V(g2)$shape <- ifelse(V(g2)$name == "A", "circle", "square")
V(g2)$color <- ifelse(V(g2)$name == "A", "lightblue", "orange")

label_expr <- expression("A", B[1], B[2])
label_expr <- label_expr[match(V(g2)$name, c("A", "B_1", "B_2"))]

plot(
  g2,
  layout = layout_v2,
  main = "Bipartite",
  vertex.size = 50,
  vertex.label = label_expr,
  vertex.label.cex = 1.5,
  vertex.frame.color = "black"
)

dev.off()
