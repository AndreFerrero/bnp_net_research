library(igraph)

set.seed(123)

# Create the graph data
edges <- data.frame(
  from = c("A", "A", "A"),
  to   = c("B", "C", "D")
)
g1 <- graph_from_data_frame(edges, directed = FALSE)

# Create a second graph with permuted labels
original_labels <- V(g1)$name
permuted_labels <- c("C", "B", "A", "D") # A and C are swapped
name_map <- setNames(permuted_labels, original_labels)

g2 <- g1
V(g2)$name <- name_map[V(g1)$name]

# Use a consistent layout for both graphs
layout <- layout_with_fr(g1)

# Function to assign consistent node colors
assign_node_colors <- function(names) {
  sapply(names, function(n) {
    if (n == "A") {
      "lightblue"
    } else if (n == "C") {
      "orange"
    } else {
      "gray90"
    }
  })
}

colors_g1 <- assign_node_colors(V(g1)$name)
colors_g2 <- assign_node_colors(V(g2)$name)

# Create directory if needed
if (!dir.exists("examples")) {
  dir.create("examples")
}

# Plot original graph (g1)
pdf("examples/node_exch_g1.pdf", width = 4, height = 2.5)
par(mar = c(0, 0, 0, 0))
plot(g1,
  layout = layout,
  vertex.label = V(g1)$name,
  vertex.label.color = "black",
  vertex.color = colors_g1,
  vertex.size = 30,
  edge.color = "gray50"
)
dev.off()

# Plot permuted graph (g2)
pdf("examples/node_exch_g2.pdf", width = 4, height = 2.5)
par(mar = c(0, 0, 0, 0))
plot(g2,
  layout = layout,
  vertex.label = V(g2)$name,
  vertex.label.color = "black",
  vertex.color = colors_g2,
  vertex.size = 30,
  edge.color = "gray50"
)
dev.off()
