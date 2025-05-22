library(readr)
library(here)
snap_folder = here("SNAP")

data_path = here(snap_folder, "gene_desease.gz")

# Read compressed TSV directly
edges <- read_tsv(data_path, col_names = T)

# assume `edges` is your 15M×2 data.frame or tibble
set.seed(123)
subn <- 1e2
edges_iid <- edges[sample(nrow(edges), size = subn), ]
# now fit your model on edges_shuffled

colnames(edges_iid) <- c("Disease", "Gene")


K_D <- unique(edges_iid$Disease) |> length()
K_G <- unique(edges_iid$Gene) |> length()

n_D <- table(edges_iid$Disease)
n_G <- table(edges_iid$Gene)

e <- nrow(edges_iid)

d_obs <- e/(K_D*K_G)

rm(edges)
rm(edges_iid)

