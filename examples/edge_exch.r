library(igraph)
set.seed(123)

# Use the same edge list as in the node exchangeability example
edges <- data.frame(
  from = c("A", "A", "A"),
  to   = c("B", "C", "D")
)

# Create original graph
g1 <- graph_from_data_frame(edges, directed = FALSE)

# Define edge colors and permutation
edge_colors <- c("red", "blue", "green")
edge_labels <- c("1", "2", "3")
edge_labels_perm <- c("3", "2", "1")

# Use same layout for both graphs
layout <- layout_with_fr(g1)
rownames(layout) <- V(g1)$name

# Compute edge label positions
compute_edge_label_positions <- function(graph, layout, d = 0.15) {
  ELx <- numeric(ecount(graph))
  ELy <- numeric(ecount(graph))

  for (i in 1:ecount(graph)) {
    ends_idx <- ends(graph, i)
    x1 <- layout[ends_idx[1], 1]
    y1 <- layout[ends_idx[1], 2]
    x2 <- layout[ends_idx[2], 1]
    y2 <- layout[ends_idx[2], 2]

    # Midpoint
    mx <- (x1 + x2) / 2
    my <- (y1 + y2) / 2

    # Compute perpendicular shift
    dx <- x2 - x1
    dy <- y2 - y1

    if (abs(dy) < 1e-3) {
      shift_x <- 0
      shift_y <- d
    } else {
      S <- -dx / dy
      shift_x <- d / sqrt(1 + S^2)
      shift_y <- S * shift_x
    }

    ELx[i] <- mx + shift_x
    ELy[i] <- my + shift_y
  }

  return(list(x = ELx, y = ELy))
}

pos_g1 <- compute_edge_label_positions(g1, layout)

# Expand plot limits
expand_limits <- function(layout, ELx, ELy, margin = 0.1) {
  x_min <- min(c(layout[, 1], ELx))
  x_max <- max(c(layout[, 1], ELx))
  y_min <- min(c(layout[, 2], ELy))
  y_max <- max(c(layout[, 2], ELy))

  xlim <- c(x_min - margin, x_max + margin)
  ylim <- c(y_min - margin, y_max + margin)

  return(list(xlim = xlim, ylim = ylim))
}

lims_g1 <- expand_limits(layout, pos_g1$x, pos_g1$y)

# Define small margins
small_margins <- c(0.5, 0, 0.5, 0)

# Plot original edge-labeled graph
pdf("examples/edge_exch_g1.pdf", width = 4, height = 2.5)
par(mar = small_margins)
plot(g1,
  layout = layout,
  vertex.label = V(g1)$name,
  vertex.label.color = "black",
  vertex.color = "gray90",
  vertex.size = 30,
  edge.label = edge_labels,
  edge.label.cex = 1.2,
  edge.color = edge_colors,
  edge.label.color = "black",
  rescale = FALSE,
  xlim = lims_g1$xlim,
  ylim = lims_g1$ylim,
  edge.label.x = pos_g1$x,
  edge.label.y = pos_g1$y,
  main = NA
)
dev.off()

# Plot permuted edge-labeled graph
pdf("examples/edge_exch_g2.pdf", width = 4, height = 2.5)
par(mar = small_margins)
plot(g1,
  layout = layout,
  vertex.label = V(g1)$name,
  vertex.label.color = "black",
  vertex.color = "gray90",
  vertex.size = 30,
  edge.label = edge_labels_perm,
  edge.label.cex = 1.2,
  edge.color = edge_colors,
  edge.label.color = "black",
  rescale = FALSE,
  xlim = lims_g1$xlim,
  ylim = lims_g1$ylim,
  edge.label.x = pos_g1$x,
  edge.label.y = pos_g1$y,
  main = NA
)
dev.off()
