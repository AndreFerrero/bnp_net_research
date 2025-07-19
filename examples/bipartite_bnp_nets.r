set.seed(123)

# Bipartite edges: pairs of nodes (A node, B node)
edges <- list(
  c(1, 1),
  c(2, 2),
  c(3, 2),
  c(2, 3),
  c(1, 1)
)

n_edges <- length(edges)

# Unique nodes in each group
nodes_A <- sort(unique(sapply(edges, `[`, 1)))
nodes_B <- sort(unique(sapply(edges, `[`, 2)))

# Colors for nodes (soft palettes)
library(RColorBrewer)
colors_A <- brewer.pal(min(length(nodes_A), 8), "Pastel1")
colors_B <- brewer.pal(min(length(nodes_B), 8), "Pastel2")
names(colors_A) <- nodes_A
names(colors_B) <- nodes_B

pdf("examples/bipartite_bnp_nets.pdf", width = 3, height = 2.5)
par(mar = c(1.5, 0, 0, 0), mgp = c(3, 0.3, 0))  # Bring labels closer

# Plot canvas
plot(NA,
  xlim = c(0.5, n_edges + 0.5), ylim = c(0, 3),
  xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", type = "n"
)

# Plot params
node_cex <- 5
get_node_radius <- function(cex) strheight("O", cex = cex)
node_radius <- get_node_radius(node_cex) / 0.71
margin <- 0.01

# X positions for edges
x_pos <- 1:n_edges
y_A <- 2
y_B <- 0.5

# Draw nodes and edges
for (i in seq_len(n_edges)) {
  edge <- edges[[i]]
  nodeA <- edge[1]
  nodeB <- edge[2]
  x <- x_pos[i]

  # A node (circle)
  points(x, y_A, pch = 21, bg = colors_A[as.character(nodeA)], cex = node_cex, lwd = 2)
  text(x, y_A, labels = nodeA, cex = 1.2, col = "black", font = 2)

  # B node (square)
  points(x, y_B, pch = 22, bg = colors_B[as.character(nodeB)], cex = node_cex, lwd = 2)
  text(x, y_B, labels = nodeB, cex = 1.2, col = "black", font = 2)

  # Edge line
  segments(
    x0 = x, y0 = y_A - node_radius + margin,
    x1 = x, y1 = y_B + node_radius - margin,
    lwd = 2, col = "black"
  )
}

# Create math expression labels: Y[1], Y[2], ...
edge_labels <- do.call(expression, lapply(1:n_edges, function(i) bquote(Y[.(i)])))

# Add axis (invisible line) with math labels close to plot
axis(1, at = x_pos, labels = edge_labels, lwd = 0)

dev.off()
