# Load required package
library(igraph)

# Updated edges for unipartite graph (including B–C)
edges1 <- c("A", "B", "A", "B", "A", "C", "B", "C")
g1 <- graph(edges1, directed = FALSE)

edges11 <- c("A", "B", "A", "B", "A", "C")
g11 <- graph(edges11, directed = FALSE)

# Bipartite graph: relabel B, C to B_1 and B_2
edges2 <- c("A", "B_1", "A", "B_1", "A", "B_2")
g2 <- graph(edges2, directed = FALSE)

g1_bin <- simplify(g1, remove.multiple = TRUE, remove.loops = TRUE)
g11_bin <- simplify(g11, remove.multiple = TRUE, remove.loops = TRUE)

# Curved edges
E(g1)$curved <- curve_multiple(g1)
E(g11)$curved <- curve_multiple(g11)
E(g2)$curved <- curve_multiple(g2)

# Define layout (V shape, A at bottom)
layout_v1 <- matrix(c(0, -1, -1, 1, 1, 1), ncol = 2, byrow = TRUE) # A, B, C
layout_v1 <- layout_v1[match(V(g1)$name, c("A", "B", "C")), ]

layout_v11 <- matrix(c(0, -1, -1, 1, 1, 1), ncol = 2, byrow = TRUE) # A, B, C
layout_v11 <- layout_v11[match(V(g11)$name, c("A", "B", "C")), ]

layout_v2 <- matrix(c(0, -1, -1, 1, 1, 1), ncol = 2, byrow = TRUE) # A, B_1, B_2
layout_v2 <- layout_v2[match(V(g2)$name, c("A", "B_1", "B_2")), ]

# Open PDF
pdf("examples/uni_net_example.pdf", width = 5, height = 5)
par(mar = c(0.5, 0.5, 2, 0.5))

# Plot 1: Unipartite with B–C edge
plot(
  g1,
  layout = layout_v1,
  vertex.size = 50,
  vertex.label.cex = 1.5,
  vertex.color = "lightblue",
  vertex.frame.color = "black"
)

dev.off()

pdf("examples/uni_net_example2_dark.pdf", width = 5, height = 5)
par(mar = c(0.5, 0.5, 2, 0.5))

# Plot 1: Unipartite with B–C edge
plot(
  g11,
  layout = layout_v11,
  vertex.size = 50,
  vertex.label.cex = 1.5,
  vertex.color = "lightblue",
  vertex.frame.color = "black",
  edge.color = "black"
)

dev.off()

pdf("examples/bip_net_example.pdf", width = 5, height = 5)
par(mar = c(0.5, 0.5, 2, 0.5))

# Plot 2: Bipartite with B[1], B[2] and shapes
V(g2)$shape <- ifelse(V(g2)$name == "A", "circle", "square")
V(g2)$color <- ifelse(V(g2)$name == "A", "orange", "orange")

label_expr <- expression("A", B[1], B[2])
label_expr <- label_expr[match(V(g2)$name, c("A", "B_1", "B_2"))]

plot(
  g2,
  layout = layout_v2,
  vertex.size = 50,
  vertex.label = label_expr,
  vertex.label.cex = 1.5,
  vertex.frame.color = "black"
)

dev.off()

pdf("examples/bin_net_example.pdf", width = 5, height = 5)
par(mar = c(0.5, 0.5, 2, 0.5))

plot(
  g1_bin,
  layout = layout_v1,
  vertex.size = 50,
  vertex.label.cex = 1.5,
  vertex.color = "lightblue",
  vertex.frame.color = "black"
)
dev.off()

pdf("examples/uni_net_example2_dark.pdf", width = 5, height = 5)

par(mar = c(0.5, 0.5, 2, 0.5))

plot(
  g11,
  layout = layout_v11,
  vertex.size = 50,
  vertex.label.cex = 1.5,
  vertex.color = "lightblue",
  vertex.frame.color = "black",
  edge.color = "darkgray"
)
dev.off()

pdf("examples/bin_net_example2_dark.pdf", width = 5, height = 5)
par(mar = c(0.5, 0.5, 2, 0.5))

plot(
  g11_bin,
  layout = layout_v11,
  vertex.size = 50,
  vertex.label.cex = 1.5,
  vertex.color = "lightblue",
  vertex.frame.color = "black",
  edge.color = "black"
)

dev.off()