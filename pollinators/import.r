# Load necessary library
library(tidyverse)
library(here)

poll_dir <- here("pollinators/mpl06005")
# Read the CSV file
# (Change the file path to where you've saved your file)
data <- read.csv(here(poll_dir, "M_PL_060_05.csv"),
  check.names = FALSE, row.names = 1
)

# Convert the adj matrix to an edge list
edges <- data |>
  rownames_to_column(var = "plant") |>
  pivot_longer(
    cols = -plant,
    names_to = "pollinator",
    values_to = "weight"
  ) |>
  filter(weight > 0) # Optional: remove zero-weight edges

plant <- edges |>
  group_by(plant) |>
  summarise(counts = sum(weight))

poll <- edges |>
  group_by(pollinator) |>
  summarise(counts = sum(weight))

(d_obs <- nrow(edges) / (nrow(plant) * nrow(poll)))

subn <- sum(edges$weight) / 2 |> round()

# 1) explode into tickets
tickets <- rep(seq_len(nrow(edges)), times = edges$weight)

# 2) sample tickets *without* replacement
sel <- sample(tickets, size = subn, replace = FALSE)

# 3) count how many times each edge got picked
sub_weights_df <- as.data.frame(table(sel), stringsAsFactors = FALSE) |>
  transmute(
    row_index = as.integer(sel),
    weight    = as.integer(Freq)
  )

# 4) join back to your edge list
sub_edges <- edges |>
  slice(sub_weights_df$row_index) |>
  mutate(weight = sub_weights_df$weight)

# View
print(sub_edges)


sub_plant <- sub_edges |>
  group_by(plant) |>
  summarise(counts = sum(weight))

sub_poll <- sub_edges |>
  group_by(pollinator) |>
  summarise(counts = sum(weight))
