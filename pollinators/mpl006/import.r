# Load necessary library
library(tidyverse)
library(here)

poll_dir <- here("pollinators/mpl006")
# Read the CSV file
# (Change the file path to where you've saved your file)
data <- read.csv(here(poll_dir, "M_PL_006.csv"),
                 check.names = FALSE, row.names = 1)

# Convert the wide matrix to a long edge list
edges <- data %>%
  rownames_to_column(var = "plant") %>%
  pivot_longer(
    cols = -plant,
    names_to = "pollinator",
    values_to = "weight"
  ) %>%
  filter(weight > 0)  # Optional: remove zero-weight edges

plant <- edges %>%
  group_by(plant) %>%
  summarise(counts = sum(weight))

poll <- edges %>%
  group_by(pollinator) %>%
  summarise(counts = sum(weight))


subn <- 100

sub_edges <- edges[sample(nrow(edges), subn, replace = T, prob = edges$weight),]
sub_weighted_edges <- sub_edges %>%
  group_by(plant, pollinator) %>%
  summarise(weight = n(), .groups = "drop") %>%
  arrange(desc(weight))


sub_plant <- sub_weighted_edges %>%
  group_by(plant) %>%
  summarise(counts = sum(weight))

sub_poll <- sub_weighted_edges %>%
  group_by(pollinator) %>%
  summarise(counts = sum(weight))
