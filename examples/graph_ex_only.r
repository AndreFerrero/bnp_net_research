library(igraph)

set.seed(123)
# --- Setup (same as before) ---
edges <- data.frame(
  from = c("A", "A", "B", "A"),
  to = c("B", "C", "C", "B"),
  label = c("1", "2", "3", "4")
)
g <- graph_from_data_frame(edges, directed = FALSE)
layout <- layout_with_fr(g) * 0.5
rownames(layout) <- c("A", "B", "C")
E(g)$curved <- curve_multiple(g) * 0.2

ELx <- rep(0, ecount(g))
ELy <- rep(0, ecount(g))
for (i in 1:ecount(g)) {
  ELx[i] <- (layout[ends(g, i)[1], 1] + layout[ends(g, i)[2], 1]) / 2
  ELy[i] <- (layout[ends(g, i)[1], 2] + layout[ends(g, i)[2], 2]) / 2
}

## Adjust perpendicular to line
d <- 0.04
for (i in 1:ecount(g)) {
  if (abs(layout[ends(g, i)[1], 2] - layout[ends(g, i)[2], 2]) < 0.1) {
    ## This avoids problems with horizontal edges
    ELy[i] <- ELy[i] + shift
  } else {
    S <- (layout[ends(g, i)[2], 1] - layout[ends(g, i)[1], 1]) /
      (layout[ends(g, i)[1], 2] - layout[ends(g, i)[2], 2])
    shift <- d / sqrt(1 + S^2)
    ELx[i] <- ELx[i] + shift
    ELy[i] <- ELy[i] + S * shift
  }
}

ELx[1] <- ELx[1]
ELx[4] <- ELx[4] * (1 - 0.05)
ELy[1] <- ELy[1] * (1 + 0.1)
ELy[4] <- ELy[4] * (1 - 0.2)

pdf("graph_ex.pdf", width = 3, height = 2.5)
par(mar = c(0, 0, 0, 0))

margin <- 0.1 # tweak this value if needed
xlim <- range(layout[, 1]) + c(-margin, margin)
ylim <- range(layout[, 2]) + c(-margin, margin)

plot(g,
  layout = layout,
  vertex.label = V(g)$name,
  vertex.label.color = "black",
  vertex.color = "gray90",
  vertex.label.cex = 0.8,
  vertex.size = 8,
  edge.color = "gray50",
  edge.label = E(g)$label,
  edge.label.color = "black",
  edge.label.cex = 1,
  rescale = FALSE,
  xlim = xlim, ylim = ylim,
  edge.label.x = ELx, edge.label.y = ELy
)

dev.off()
