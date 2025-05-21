library(igraph)

# 1. Build a bipartite graph from the edges data frame
nodes <- unique(c(edges$Disease, edges$Gene))
node_types <- ifelse(nodes %in% edges$Disease, TRUE, FALSE)

g <- graph_from_data_frame(edges, vertices = data.frame(name = nodes, type = node_types), directed = FALSE)

# 2. Basic bipartite layout visualization
V(g)$color <- ifelse(V(g)$type, "skyblue", "salmon")  # Different color for each type
V(g)$size <- 5
V(g)$label <- NA  # Optional: remove labels for clarity

plot(g, layout = layout_as_bipartite(g),
     vertex.label.cex = 0.7,
     main = "Disease–Gene Bipartite Network")


# SUBSET

library(igraph)

# Calculate disease degrees (number of connected genes)
disease_degrees <- table(edges$Disease)
top_diseases <- names(sort(disease_degrees, decreasing = TRUE))[1:10]

# Subset the edges to only include top diseases
sub_edges <- edges[edges$Disease %in% top_diseases, ]

# Create graph
nodes <- unique(c(sub_edges$Disease, sub_edges$Gene))
node_types <- ifelse(nodes %in% sub_edges$Disease, TRUE, FALSE)

g_sub <- graph_from_data_frame(sub_edges, vertices = data.frame(name = nodes, type = node_types), directed = FALSE)

# Plot with bipartite layout
V(g_sub)$color <- ifelse(V(g_sub)$type, "skyblue", "salmon")
V(g_sub)$size <- 5
V(g_sub)$label <- NA

plot(g_sub, layout = layout_as_bipartite, main = "Top 50 Disease–Gene Bipartite Network")
