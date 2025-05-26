library(readr)
library(here)
micro_folder <- here("micro")

data_path <- here(micro_folder, "HMP_132902142.csv")
edges <- read_csv(data_path, col_select = c(1,2))

# assume `edges` is your 15M×2 data.frame or tibble
set.seed(123)
subn <- 500
edges_iid <- edges[sample(nrow(edges), size = subn), ]
# now fit your model on edges_shuffled

colnames(edges_iid) <- c("bacteria", "location")


K_A <- unique(edges_iid$bacteria) |> length()
K_B <- unique(edges_iid$location) |> length()

n_A <- table(edges_iid$bacteria)
n_B <- table(edges_iid$location)

e <- nrow(edges_iid)

d_obs <- e/(K_A*K_B)

rm(edges)
rm(edges_iid)