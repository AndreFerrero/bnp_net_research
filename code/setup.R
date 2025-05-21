
library(igraph)
library(bmotif)
library(ggplot2)
library(viridis)
library(dplyr)
library(here)
library(qgraph)
library(parallel)

functions_folder = here("code/functions")
funs <- list.files(functions_folder, pattern = "\\.R$", full.names = TRUE)

for(f in funs){
  source(f)
}

rm(f)
rm(funs)


