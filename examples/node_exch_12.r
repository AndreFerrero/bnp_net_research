library(igraph)

set.seed(123)

edges <- data.frame(
  from = c("A", "A", "A"),
  to   = c("B", "C", "D")
)
g1 <- graph_from_data_frame(edges, directed = FALSE)

original_labels <- V(g1)$name
permuted_labels <- c("C", "B", "A", "D")
name_map <- setNames(permuted_labels, original_labels)

g2 <- g1
V(g2)$name <- name_map[V(g1)$name]

layout <- layout_with_fr(g1)

changed_labels <- V(g1)$name != V(g2)$name
label_colors_g2 <- ifelse(changed_labels, "red", "black")

# Smaller margins around plot
small_margins <- c(1, 0, 1, 0) # bottom, left, top, right

# Plot original graph without title and small margins
pdf("examples/node_exch_g1.pdf", width = 4, height = 2.5)
par(mar = small_margins)
plot(g1,
  layout = layout,
  vertex.label = V(g1)$name,
  vertex.label.color = "black",
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  main = NA
)
dev.off()

# Plot permuted graph without title and small margins
pdf("examples/node_exch_g2.pdf", width = 4, height = 2.5)
par(mar = small_margins)
plot(g2,
  layout = layout,
  vertex.label = V(g2)$name,
  vertex.label.color = label_colors_g2,
  vertex.color = "gray90",
  vertex.size = 30,
  edge.color = "gray50",
  main = NA
)
dev.off()
