library(igraph)
set.seed(123)

# --- Define simplified edge list with 4 edges ---
edge_sequence <- data.frame(
  from = c("A", "A", "B", "C"),
  to   = c("B", "C", "C", "D")
)

# Create original graph
g1 <- graph_from_data_frame(edge_sequence, directed = FALSE)

# Define a color for each edge
edge_colors <- c("red", "blue", "green", "orange")

# --- Permute the edge order: move first edge (A,B) to the end ---
new_order <- c(2:4, 1) # Move (A,B) to the end
edge_sequence_perm <- edge_sequence[new_order, ]

# Create permuted graph
g2 <- graph_from_data_frame(edge_sequence_perm, directed = FALSE)

# --- Use same layout for consistency ---
layout <- layout_with_fr(g1)
rownames(layout) <- c("A", "B", "C", "D")
# --- Functions to get edge colors and labels ---
get_edge_colors <- function(graph, edge_df, color_vec) {
  el <- as_edgelist(graph)
  colors <- character(nrow(el))
  for (i in seq_len(nrow(el))) {
    for (j in seq_len(nrow(edge_df))) {
      if (setequal(el[i, ], as.character(edge_df[j, ]))) {
        colors[i] <- color_vec[j]
        break
      }
    }
  }
  return(colors)
}

get_edge_labels <- function(graph, edge_df) {
  el <- as_edgelist(graph)
  labels <- character(nrow(el))
  for (i in seq_len(nrow(el))) {
    for (j in seq_len(nrow(edge_df))) {
      if (setequal(el[i, ], as.character(edge_df[j, ]))) {
        labels[i] <- as.character(j)
        break
      }
    }
  }
  return(labels)
}

# Get consistent edge colors and labels for both graphs
edge_colors_g1 <- get_edge_colors(g1, edge_sequence, edge_colors)
edge_colors_g2 <- get_edge_colors(g2, edge_sequence_perm, edge_colors)

edge_labels_g1 <- get_edge_labels(g1, edge_sequence)
edge_labels_g2 <- get_edge_labels(g2, edge_sequence_perm)

# --- Calculate manually shifted edge label positions ---
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
      # Near horizontal edge: shift vertically by +d
      shift_x <- 0
      shift_y <- d
    } else {
      # Slope of perpendicular line
      S <- -dx / dy
      shift_x <- d / sqrt(1 + S^2)
      shift_y <- S * shift_x
    }

    ELx[i] <- mx + shift_x
    ELy[i] <- my + shift_y
  }

  return(list(x = ELx, y = ELy))
}

# Compute positions for both graphs
pos_g1 <- compute_edge_label_positions(g1, layout)
pos_g2 <- compute_edge_label_positions(g2, layout)

# --- Dynamically expand plot limits based on layout and shifted edge labels ---
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
lims_g2 <- expand_limits(layout, pos_g2$x, pos_g2$y)

# --- Plot both graphs ---
pdf("edge_exch.pdf", width = 4, height = 2.5)
par(mfrow = c(1, 2), mar = c(0, 0, 2, 0))

plot(g1,
  layout = layout,
  vertex.label = V(g1)$name,
  vertex.label.color = "black",
  vertex.color = "gray90",
  vertex.size = 30,
  edge.label = edge_labels_g1,
  edge.label.cex = 1.2,
  edge.color = edge_colors_g1,
  edge.label.color = "black",
  rescale = FALSE,
  xlim = lims_g1$xlim,
  ylim = lims_g1$ylim,
  edge.label.x = pos_g1$x,
  edge.label.y = pos_g1$y,
  main = "(a)"
)

plot(g2,
  layout = layout,
  vertex.label = V(g2)$name,
  vertex.label.color = "black",
  vertex.color = "gray90",
  vertex.size = 30,
  edge.label = edge_labels_g2,
  edge.label.cex = 1.2,
  edge.color = edge_colors_g2,
  edge.label.color = "black",
  rescale = FALSE,
  xlim = lims_g2$xlim,
  ylim = lims_g2$ylim,
  edge.label.x = pos_g2$x,
  edge.label.y = pos_g2$y,
  main = "(b)"
)

dev.off()
