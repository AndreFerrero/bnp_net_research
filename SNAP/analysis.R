# =============================
# Bipartite Network Analysis: Disease–Gene Network
# =============================

# ---- Load Required Libraries ----
library(readr)      # For fast file reading
library(igraph)     # For graph construction and network metrics
library(Matrix)     # For handling sparse matrices if needed
library(bipartite)  # For ecological/bipartite network metrics

# ---- Step 1: Import the Data ----
# Change the file path as appropriate.
# Here we read the first 5,000 rows for demonstration.
file_path <- "C:/Users/anferrer/OneDrive - Alma Mater Studiorum Università di Bologna/University/Internship/research/SNAP/gene_desease.gz"
edges <- read_tsv(file_path, col_names = FALSE, n_max = 5000)
colnames(edges) <- c("Disease", "Gene")
cat("Imported", nrow(edges), "edges.\n")

# ---- Step 2: Create the Bipartite Graph ----
# All nodes (diseases and genes)
all_nodes <- unique(c(edges$Disease, edges$Gene))
# Mark type: TRUE for disease, FALSE for gene.
node_types <- ifelse(all_nodes %in% edges$Disease, TRUE, FALSE)

# Construct the bipartite graph
g <- graph_from_data_frame(edges,
                           vertices = data.frame(name = all_nodes, type = node_types),
                           directed = FALSE)
cat("Graph has", vcount(g), "nodes and", ecount(g), "edges.\n")

# ---- Step 3: Basic Bipartite Metrics ----
disease_nodes <- V(g)[V(g)$type == TRUE]
gene_nodes    <- V(g)[V(g)$type == FALSE]
max_possible_edges <- length(disease_nodes) * length(gene_nodes)
connectance <- ecount(g) / max_possible_edges
cat("Bipartite connectance (density):", round(connectance, 5), "\n")

# Incidence matrix creation
inc_mat <- table(edges$Disease, edges$Gene)
net_level <- networklevel(inc_mat)
cat("Network-level metrics (from bipartite::networklevel()):\n")
print(net_level)

# ---- Step 4: Degree Distributions ----
deg_disease <- igraph::degree(g, v = disease_nodes)
deg_gene    <- igraph::degree(g, v = gene_nodes)
cat("Summary of Disease Degrees:\n")
print(summary(deg_disease))
cat("Summary of Gene Degrees:\n")
print(summary(deg_gene))

# Plot degree distributions
par(mfrow = c(1,2))
hist(deg_disease, breaks = 100, main = "Disease Degree Distribution", xlab = "Degree", col = "lightblue")
hist(deg_gene, breaks = 100, main = "Gene Degree Distribution", xlab = "Degree", col = "salmon")
par(mfrow = c(1,1))

# ---- Step 4.1: Cross-Type Degree Assortativity ----
cat("Computing cross-type degree assortativity...\n")

deg_all <- igraph::degree(g)  # <-- ensure correct function is used
edge_list <- as_edgelist(g)
type_vec <- V(g)$type  # TRUE = disease, FALSE = gene

gene_degs <- ifelse(type_vec[edge_list[, 1]], deg_all[edge_list[, 2]], deg_all[edge_list[, 1]])
disease_degs <- ifelse(!type_vec[edge_list[, 1]], deg_all[edge_list[, 2]], deg_all[edge_list[, 1]])

cross_assortativity <- cor(gene_degs, disease_degs, method = "pearson")
cat("Cross-type degree assortativity (gene–disease correlation):", round(cross_assortativity, 4), "\n")


# ---- Step 5: Centrality Measures ----
cat("Calculating betweenness centrality (this may take time)...\n")
betw_disease <- betweenness(g, v = disease_nodes)
betw_gene    <- betweenness(g, v = gene_nodes)
cat("Top 5 diseases by betweenness:\n")
print(head(sort(betw_disease, decreasing = TRUE), 5))
cat("Top 5 genes by betweenness:\n")
print(head(sort(betw_gene, decreasing = TRUE), 5))

clo_disease <- closeness(g, vids = disease_nodes)
clo_gene    <- closeness(g, vids = gene_nodes)
cat("Summary closeness for diseases:\n")
print(summary(clo_disease))
cat("Summary closeness for genes:\n")
print(summary(clo_gene))

# ---- Step 6: Component Analysis ----
comp <- igraph::components(g)
cat("Number of components in the bipartite network:", comp$no, "\n")
cat("Component size distribution:\n")
print(table(comp$csize))

giant_comp_id <- which.max(comp$csize)
g_giant <- induced_subgraph(g, which(comp$membership == giant_comp_id))
cat("Giant component contains", vcount(g_giant), "nodes and", ecount(g_giant), "edges.\n")

# ---- Step 7: Projections to One-Mode Networks ----
proj <- bipartite_projection(g_giant)
g_disease <- simplify(proj$proj1, remove.loops = TRUE)
g_gene    <- simplify(proj$proj2, remove.loops = TRUE)

cat("Disease projection: ", vcount(g_disease), "nodes,", ecount(g_disease), "edges.\n")
cat("Gene projection: ", vcount(g_gene), "nodes,", ecount(g_gene), "edges.\n")

# Plot degree distribution in the gene projection
deg_gene_proj <- igraph::degree(g_gene)
hist(deg_gene_proj, breaks = 100, main = "Gene–Gene Projection Degree Distribution", xlab = "Degree", col = "lightgreen")

# ---- Step 8: Community Detection on the Gene Projection ----
cat("Running community detection (Louvain) on gene projection...\n")
comm_gene <- cluster_louvain(g_gene)
cat("Number of communities detected:", length(comm_gene), "\n")
cat("Sizes of communities:\n")
print(sizes(comm_gene))
barplot(sizes(comm_gene), main = "Community Sizes in Gene–Gene Projection", ylab = "Number of Genes", col = "gray80")

V(g_gene)$community <- membership(comm_gene)

# ---- Step 9: Additional Bipartite Metrics (Nestedness and Specialization) ----
nested_val <- nested(inc_mat, method = "NODF")
cat("Nestedness (NODF):", round(nested_val, 4), "\n")

specialization <- networklevel(inc_mat, index = "H2")
cat("Network-level specialization (H2):", round(specialization$H2, 4), "\n")

# ---- Final Summary Output ----
cat("\n=== Analysis Summary ===\n")
cat("Bipartite Connectance:", round(connectance, 5), "\n")
cat("Network-level metrics from bipartite::networklevel():\n")
print(net_level)
cat("Nestedness (NODF):", round(nested_val, 4), "\n")
cat("Specialization (H2):", round(specialization$H2, 4), "\n")
cat("Number of components:", comp$no, "\n")
cat("Giant component size:", vcount(g_giant), "nodes\n")
cat("Gene projection: ", vcount(g_gene), "nodes and", ecount(g_gene), "edges.\n")
cat("Detected communities in gene projection:", length(comm_gene), "\n")
cat("Cross-type degree assortativity:", round(cross_assortativity, 4), "\n")
cat("\nAnalysis complete.\n")
