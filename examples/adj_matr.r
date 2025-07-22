library(igraph)

# Graphs with edges (including B-C in unipartite)
edges1 <- c("A", "B", "A", "B", "A", "C", "B", "C")
g1 <- graph(edges1, directed = FALSE)

edges2 <- c("A", "B_1", "A", "B_1", "A", "B_2")
g2 <- graph(edges2, directed = FALSE)

adj1 <- as.matrix(as_adjacency_matrix(g1, sparse = FALSE))
adj2 <- as.matrix(as_adjacency_matrix(g2, sparse = FALSE))

adj1 <- adj1[c("A", "B", "C"), c("A", "B", "C")]
adj2 <- adj2[c("A", "B_1", "B_2"), c("A", "B_1", "B_2")]

# Custom colors: white (no edge) to blue (edge)
col_palette_blue <- colorRampPalette(c("white", "dodgerblue4"))(2)
col_palette_orange <- colorRampPalette(c("white", "darkorange3"))(2)


# Function to plot adjacency matrix nicely with optional expression labels

plot_adj_matrix <- function(adj, use_expression_labels = FALSE, palette = col_palette_blue) {
  n <- nrow(adj)

  image(
    1:n, 1:n, t(adj[n:1, ]),
    col = palette,
    axes = FALSE,
    xlab = "", ylab = "",
    xlim = c(0.5, n + 0.5),
    ylim = c(0.5, n + 0.5),
    asp = 1
  )

  abline(h = 0.5:(n + 0.5), v = 0.5:(n + 0.5), col = "grey80")

  if (use_expression_labels) {
    x_labels <- c(expression(A[phantom(1)]), expression(B[1]), expression(B[2]))
    y_labels <- rev(x_labels)
    axis(1, at = 1:n, labels = x_labels, tick = FALSE, cex.axis = 1.2)
    axis(2, at = 1:n, labels = y_labels, tick = FALSE, las = 2, cex.axis = 1.2)
  } else {
    axis(1, at = 1:n, labels = colnames(adj), tick = FALSE, cex.axis = 1.2)
    axis(2, at = 1:n, labels = rev(rownames(adj)), tick = FALSE, las = 2, cex.axis = 1.2)
  }

  for (i in 1:n) {
    for (j in 1:n) {
      val <- adj[i, j]
      if (val > 0) {
        text(j, n - i + 1, labels = val, col = "black", cex = 1.3, font = 2)
      }
    }
  }
}


# Unipartite: blue
pdf("examples/adj_matr_uni.pdf", width = 3.5, height = 3.5)
par(mar = c(3, 3, 1, 1))
plot_adj_matrix(adj1, use_expression_labels = FALSE, palette = col_palette_blue)
dev.off()

# Bipartite: orange
pdf("examples/adj_matr_bip.pdf", width = 3.5, height = 3.5)
par(mar = c(3, 3, 1, 1))
plot_adj_matrix(adj2, use_expression_labels = TRUE, palette = col_palette_orange)
dev.off()
