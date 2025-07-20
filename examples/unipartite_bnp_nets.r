set.seed(123)

# Define edges: pairs of node labels (integers)
# Some labels repeat to show identical values from DP
edges <- list(
  c(1, 2),
  c(3, 4),
  c(2, 5),
  c(3, 6),
  c(2, 5)
)

n_edges <- length(edges)

# Extract all node labels
all_nodes <- unique(unlist(edges))

# Assign colors (soft palette)
library(RColorBrewer)
node_colors <- brewer.pal(min(length(all_nodes), 8), "Pastel1")
names(node_colors) <- all_nodes

pdf("examples/unipartite_bnp_nets.pdf", width = 3, height = 2.5)
par(mar = c(1.5, 2, 0, 0.5), mgp = c(3, 0.3, 0)) # reduced bottom margin, no left/top/right margins

plot(NA,
  xlim = c(0.5, n_edges + 0.5), ylim = c(0, 3),
  xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", type = "n"
)

# Plot params
node_cex <- 5
get_node_radius <- function(cex) strheight("O", cex = cex)
node_radius <- get_node_radius(node_cex)/0.72
margin <- 0.01

# X positions for edges
x_pos <- 1:n_edges
y_top <- 2
y_bottom <- 0.5

# Add "G:" label on the left, centered vertically
text(x = 0, y = mean(c(y_top, y_bottom)), labels = expression(G), font = 2, cex = 1.4, xpd = NA)

# Draw edges and nodes
for (i in seq_len(n_edges)) {
  nodes_pair <- edges[[i]]
  x <- x_pos[i]

  # Draw top node (circle)
  points(x, y_top, pch = 21, bg = node_colors[as.character(nodes_pair[1])], cex = node_cex, lwd = 2)
  text(x, y_top, labels = nodes_pair[1], cex = 1.2, col = "black", font = 2)

  # Draw bottom node (circle)
  points(x, y_bottom, pch = 21, bg = node_colors[as.character(nodes_pair[2])], cex = node_cex, lwd = 2)
  text(x, y_bottom, labels = nodes_pair[2], cex = 1.2, col = "black", font = 2)

  # Draw edge line
  segments(
    x0 = x, y0 = y_top - node_radius + margin,
    x1 = x, y1 = y_bottom + node_radius - margin,
    lwd = 2, col = "black"
  )
}

# Create expression labels for the axis: Y_1, Y_2, ..., Y_n
edge_labels <- do.call(expression, lapply(1:n_edges, function(i) bquote(Y[.(i)])))

# Add x-axis without axis line (lwd=0), with expression labels
axis(1, at = x_pos, labels = edge_labels, lwd = 0, cex.axis = 1.4)

dev.off()
